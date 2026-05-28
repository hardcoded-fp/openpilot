#!/usr/bin/env python3

import ast
import datetime
import html
import logging
import os
import pprint
import shlex
import subprocess
import tempfile
import time

logging.basicConfig(level=logging.INFO, format="%(message)s")

OPENPILOT_REPO_URL = "https://github.com/commaai/openpilot.git"
SUPPORTED_BASE_BRANCHES = ("nightly", "nightly-dev", "release-mici", "release-tizi")
UPSTREAM_BRANCH_CANDIDATES = {
    "nightly": ("nightly",),
    "nightly-dev": ("nightly-dev",),
    "release-mici": ("release-mici",),
    "release-tizi": ("release-tizi",),
}
TRANSIENT_GIT_ERROR_MARKERS = (
    "The requested URL returned error: 500",
    "The requested URL returned error: 502",
    "The requested URL returned error: 503",
    "The requested URL returned error: 504",
    "Internal Server Error",
    "remote end hung up unexpectedly",
    "Connection reset",
    "Connection timed out",
    "Failed to connect",
    "Operation timed out",
    "TLS connection was non-properly terminated",
    "Could not resolve host",
    "early EOF",
    "RPC failed",
)
TRANSIENT_GIT_RETURN_CODES = (
    -15,  # SIGTERM when subprocess reports the signal directly.
    143,  # 128 + SIGTERM when reported through the shell.
)
PUSH_RETRY_ATTEMPTS = 5
PUSH_RETRY_BASE_DELAY_SECONDS = 10
PUSH_BATCH_SIZE = 100
BRANCH_CONTEXTS = {}


def run(cmd, check=True):
    logging.info("+ %s", cmd)
    return subprocess.run(cmd, shell=True, check=check)


def capture(cmd):
    logging.info("+ %s", cmd)
    return subprocess.check_output(cmd, shell=True, text=True).strip()


def run_for_output(cmd, check=True):
    logging.info("+ %s", cmd)
    process = subprocess.run(
        cmd,
        shell=True,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if process.stdout:
        logging.info(process.stdout.rstrip())
    if check and process.returncode != 0:
        raise subprocess.CalledProcessError(
            process.returncode, process.args, output=process.stdout
        )
    return process


def run_args_for_output(args, cwd=None, check=True, log_output=True):
    logging.info("+ %s", " ".join(shlex.quote(str(arg)) for arg in args))
    process = subprocess.run(
        args,
        cwd=cwd,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if log_output and process.stdout:
        logging.info(process.stdout.rstrip())
    if check and process.returncode != 0:
        raise subprocess.CalledProcessError(
            process.returncode, process.args, output=process.stdout
        )
    return process


def is_transient_git_error(process):
    if process.returncode in TRANSIENT_GIT_RETURN_CODES:
        return True

    output = process.stdout or ""
    normalized_output = output.lower()
    return any(
        marker.lower() in normalized_output for marker in TRANSIENT_GIT_ERROR_MARKERS
    )


def git_cmd(args):
    return "git " + " ".join(shlex.quote(str(arg)) for arg in args)


def git_capture(args, cwd="comma_openpilot", env=None):
    logging.info("+ cd %s && %s", cwd, git_cmd(args))
    return subprocess.check_output(
        ["git", *args], cwd=cwd, env=env, text=True
    ).strip()


def git_capture_bytes(args, cwd="comma_openpilot", env=None):
    logging.info("+ cd %s && %s", cwd, git_cmd(args))
    return subprocess.check_output(["git", *args], cwd=cwd, env=env)


def git_run(args, cwd="comma_openpilot", env=None):
    logging.info("+ cd %s && %s", cwd, git_cmd(args))
    return subprocess.run(["git", *args], cwd=cwd, env=env, check=True)


def get_branch_context(upstream_branch):
    if upstream_branch in BRANCH_CONTEXTS:
        return BRANCH_CONTEXTS[upstream_branch]

    base_ref = f"origin/{upstream_branch}"
    launch_env_meta = git_capture(["ls-tree", base_ref, "launch_env.sh"])
    context = {
        "base_ref": base_ref,
        "base_commit": git_capture(["rev-parse", base_ref]),
        "base_tree": git_capture(["rev-parse", f"{base_ref}^{{tree}}"]),
        "launch_env_mode": launch_env_meta.split()[0],
        "launch_env": git_capture_bytes(["show", f"{base_ref}:launch_env.sh"]),
        "commit_date": git_capture(["log", "-1", "--format=%cd", "--date=iso-strict", base_ref]),
        "author_date": git_capture(["log", "-1", "--format=%ad", "--date=iso-strict", base_ref]),
    }
    BRANCH_CONTEXTS[upstream_branch] = context
    return context


def push_branch_with_retries(branch_name):
    push_branch_batch_with_retries([branch_name])


def push_branch_batch_with_retries(branch_names):
    push_args = ["git", "push", "origin", "--force"]
    push_args.extend(f"{branch_name}:{branch_name}" for branch_name in branch_names)
    for attempt in range(1, PUSH_RETRY_ATTEMPTS + 1):
        process = run_args_for_output(push_args, cwd="comma_openpilot", check=False)
        if process.returncode == 0:
            return

        output = process.stdout or ""
        if not is_transient_git_error(process) or attempt == PUSH_RETRY_ATTEMPTS:
            raise subprocess.CalledProcessError(
                process.returncode, process.args, output=output
            )

        delay = PUSH_RETRY_BASE_DELAY_SECONDS * attempt
        logging.warning(
            "Push batch of %s branches failed with a transient GitHub error. "
            "Retrying in %s seconds (%s/%s). Return code: %s.",
            len(branch_names),
            delay,
            attempt + 1,
            PUSH_RETRY_ATTEMPTS,
            process.returncode,
        )
        time.sleep(delay)


def chunked(items, size):
    for index in range(0, len(items), size):
        yield items[index : index + size]


def get_remote_branch_tips(branch_names):
    remote_tips = {}
    for batch in chunked(branch_names, PUSH_BATCH_SIZE):
        patterns = [f"refs/heads/{branch_name}" for branch_name in batch]
        process = run_args_for_output(
            ["git", "ls-remote", "origin", *patterns],
            cwd="comma_openpilot",
            check=False,
            log_output=False,
        )
        if process.returncode != 0:
            raise subprocess.CalledProcessError(
                process.returncode, process.args, output=process.stdout
            )

        for line in process.stdout.splitlines():
            sha, ref = line.split("\t", 1)
            remote_tips[ref.removeprefix("refs/heads/")] = sha
    return remote_tips


def get_local_branch_tips(branch_names):
    local_tips = {}
    for batch in chunked(branch_names, PUSH_BATCH_SIZE):
        refs = [f"refs/heads/{branch_name}" for branch_name in batch]
        output = git_capture(
            ["for-each-ref", "--format=%(refname:short) %(objectname)", *refs]
        )
        if not output:
            continue

        for line in output.splitlines():
            branch_name, sha = line.split(" ", 1)
            local_tips[branch_name] = sha
    return local_tips


def get_available_upstream_branches():
    refs = capture(f"git ls-remote --heads {OPENPILOT_REPO_URL}")
    branches = set()
    for line in refs.splitlines():
        _, ref = line.split("\t", 1)
        branches.add(ref.removeprefix("refs/heads/"))
    return branches


def get_string_literal(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def resolve_base_branches():
    available_branches = get_available_upstream_branches()
    resolved = {}
    skipped = []
    for public_branch in SUPPORTED_BASE_BRANCHES:
        for candidate in UPSTREAM_BRANCH_CANDIDATES[public_branch]:
            if candidate in available_branches:
                resolved[public_branch] = candidate
                break
        if public_branch not in resolved:
            skipped.append(public_branch)

    if not resolved:
        raise RuntimeError(
            "Failed to resolve any supported upstream branches. "
            + f". Available branches: {', '.join(sorted(available_branches))}"
        )

    logging.info("Resolved base branches: %s", pprint.pformat(resolved))
    if skipped:
        logging.info("Skipping unavailable branches: %s", ", ".join(skipped))
    return resolved


def parse_cars(upstream_branch):
    """
    Parse car names from the current upstream opendbc layout.

    We read the `CAR` class values from `values.py` and fall back to
    `fingerprints.py` keys for branches that still define models there.

    Example:

    class CAR:
      # Hyundai
      ELANTRA = "HYUNDAI ELANTRA 2017"
      ELANTRA_2021 = "HYUNDAI ELANTRA 2021"
      ELANTRA_HEV_2021 = "HYUNDAI ELANTRA HYBRID 2021"
      HYUNDAI_GENESIS = "HYUNDAI GENESIS 2015-2016"
      IONIQ = "HYUNDAI IONIQ HYBRID 2017-2019"

    Another format is to look in `fingerprints.py` instead of `values.py`.

    Example:

    FW_VERSIONS = {
        CAR.TOYOTA_AVALON: {
        (Ecu.abs, 0x7b0, None): [
            b'F152607060\x00\x00\x00\x00\x00\x00',
        ],
    ... continued }

    `cars` should be an array of strings like this:

    [
      "HYUNDAI ELANTRA 2017",
      "HYUNDAI ELANTRA 2021",
      "HYUNDAI ELANTRA HYBRID 2021",
      "HYUNDAI GENESIS 2015-2016",
      "HYUNDAI IONIQ HYBRID 2017-2019"
      ...
    ]

    For the second one, it should be like this:

    [
        "TOYOTA AVALON",
        ...
    ]
    """
    # Checkout branch
    run(f"cd comma_openpilot && git checkout --force origin/{upstream_branch}")

    car_root = "comma_openpilot/opendbc/car"
    values_py_paths = []
    fingerprints_py_paths = []
    cars = []

    for root, dirs, files in os.walk(car_root):
        values_py_paths += [os.path.join(root, f) for f in files if f == "values.py"]
        fingerprints_py_paths += [os.path.join(root, f) for f in files if f == "fingerprints.py"]

    for path in values_py_paths:
        logging.info("Parsing %s", path)
        with open(path, "r") as f:
            tree = ast.parse(f.read())
            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef) and node.name == "CAR":
                    for c in node.body:
                        if isinstance(c, ast.Assign):
                            string_value = get_string_literal(c.value)
                            if string_value is not None:
                                cars.append(string_value)
                            elif isinstance(c.targets[0], ast.Name):
                                # opendbc has most of the cars refactor
                                cars.append(c.targets[0].id)
                            # Sometimes it's an object initializer,
                            # If so, use the first argument
                            elif isinstance(c.value, ast.Call):
                                # Sometimes
                                if len(c.value.args) > 0:
                                    string_value = get_string_literal(c.value.args[0])
                                    if string_value is not None:
                                        cars.append(string_value)

    for path in fingerprints_py_paths:
        logging.info("Parsing %s", path)
        with open(path, "r") as f:
            tree = ast.parse(f.read())
            for node in ast.walk(tree):
                if isinstance(node, ast.Assign) and isinstance(node.value, ast.Dict):
                    for key in node.value.keys:
                        if isinstance(key, ast.Attribute):
                            cars.append(key.attr)

    cars = list(dict.fromkeys(cars))
    logging.info("Found %d cars in %s", len(cars), upstream_branch)

    return cars


def prepare_op_repo():
    """
    Prepare the openpilot repo and resolve the current upstream deployment branches.
    """
    resolved_branches = resolve_base_branches()

    # Try to clone the repo to comma_openpilot.
    # If it fails, it means it already exists, so we can just pull
    # the latest changes.
    logging.info("Setting up openpilot repo. Ignore errors if it already exists.")

    if not os.path.isdir("comma_openpilot/.git"):
        run(
            "git clone --filter=blob:none --no-checkout "
            f"{OPENPILOT_REPO_URL} comma_openpilot"
        )
    # Make sure that comma_openpilot is usiing that as the origin.
    run(f"cd comma_openpilot && git remote set-url origin {OPENPILOT_REPO_URL}")
    branch_refspecs = " ".join(
        f"+refs/heads/{branch}:refs/remotes/origin/{branch}"
        for branch in sorted(set(resolved_branches.values()))
    )
    run(
        "cd comma_openpilot && git fetch --prune --depth 1 origin "
        f"{branch_refspecs}"
    )

    logging.info("Done setting up openpilot repo.")
    return resolved_branches


def sanitize_branch_component(car):
    return (
        car.replace(" ", "-")
        .replace("&", "AND")
        .replace("(", "")
        .replace(")", "")
        .replace("_", "-")
        .lower()
    )


def generate_branch(base, upstream_branch, car):
    """
    Make a new branch for the car with a hardcoded fingerprint
    """

    # - instead of _ because the keyboard is one tap for - vs two for _
    # & is AND because & may be too special
    # Lowercase because there's no caps lock in the keyboard
    # Remove () because they are special characters and may cause issues
    branch_name = f"{base}-{sanitize_branch_component(car)}"
    logging.info("Generating branch %s from origin/%s", branch_name, upstream_branch)

    branch_context = get_branch_context(upstream_branch)
    base_commit = branch_context["base_commit"]
    base_tree = branch_context["base_tree"]
    launch_env_mode = branch_context["launch_env_mode"]
    launch_env = branch_context["launch_env"]
    launch_env += f'\nexport FINGERPRINT="{car}"\n'.encode()
    # hash-object needs stdin; keep this separate so log output stays readable.
    logging.info("+ cd comma_openpilot && git hash-object -w --stdin < launch_env.sh")
    launch_env_blob = subprocess.check_output(
        ["git", "hash-object", "-w", "--stdin"],
        cwd="comma_openpilot",
        input=launch_env,
        text=False,
    ).decode().strip()

    commit_date = branch_context["commit_date"]
    author_date = branch_context["author_date"]
    index_file = tempfile.NamedTemporaryFile(prefix="hardcoded-fp-index-", delete=False)
    index_file.close()

    try:
        env = os.environ.copy()
        env["GIT_INDEX_FILE"] = index_file.name
        git_run(["read-tree", base_tree], env=env)
        git_run(
            [
                "update-index",
                "--add",
                "--cacheinfo",
                f"{launch_env_mode},{launch_env_blob},launch_env.sh",
            ],
            env=env,
        )
        tree = git_capture(["write-tree"], env=env)

        commit_env = env.copy()
        commit_env["GIT_AUTHOR_DATE"] = author_date
        commit_env["GIT_COMMITTER_DATE"] = commit_date
        commit = git_capture(
            [
                "commit-tree",
                tree,
                "-p",
                base_commit,
                "-m",
                f"Hardcode fingerprint for {car}",
            ],
            env=commit_env,
        )
        git_run(["update-ref", f"refs/heads/{branch_name}", commit])
    finally:
        os.unlink(index_file.name)

    return branch_name


def generate_html(base_cars_base_branches):
    # Generate a date for the page
    now = datetime.datetime.now()
    now_str = now.strftime("%Y-%m-%d %H:%M:%S UTC")

    header = """
<html>
<head>
<title>Hardcoded Fingerprint comma.ai openpilot Continuous Micro-Fork Generator branches</title>
<style>
body {
    font-family: sans-serif;
    margin: 24px;
    line-height: 1.4;
}
input[type="search"] {
    width: min(100%, 420px);
    padding: 8px 10px;
    font: inherit;
}
.filter-bar {
    margin: 20px 0 24px;
}
.filter-meta {
    color: #555;
    margin-top: 8px;
}
.branch-section {
    margin-top: 24px;
}
.car-item[hidden],
.branch-section[hidden] {
    display: none;
}
</style>
<link href="data:image/x-icon;base64,AAABAAEAEBAQAAEABAAoAQAAFgAAACgAAAAQAAAAIAAAAAEABAAAAAAAgAAAAAAAAAAAAAAAEAAAAAAAAAAAAAAA6OjoAP///wAAAOMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAzMwAAAAAAAzMzMAAAAAADMzMzAAAAAAMjMzMwAAAAAyIzMzAAAAAAMiMzAiAAAAADIjAhEgAAAAAzAiERIAAAAAAiIhESAAAAACIiIREAAAAAAiIiEgAAAAAAIiIiAAAAAAACIiAAAAAAAAAAAAD//wAAw/8AAIH/AAAA/wAAAH8AAAA/AAAAHwAAgA8AAMAHAADgAwAA8AEAAPgBAAD8AQAA/gEAAP8DAAD/hwAA" rel="icon" type="image/x-icon" />
</head>
<body>
<h1>Hardcoded Fingerprint comma.ai openpilot Continuous Micro-Fork Generator branches</h1>
<p>
<em>PRESCRIPTION ONLY: Consult your vehicle brand's <a href="https://discord.comma.ai">Discord channel</a> for guidance first.</em>
</p>
<p>
⚠️ Only to be used as a last resort! ⚠️
</p>
<p>
""" + f"""
This page was generated on {now_str}.
</p>
<p>
This is a list of all the branches with hardcoded fingerprints generated by the <a href="https://github.com/hardcoded-fp/openpilot/"> Hardcoded Fingerprint comma.ai openpilot Continuous Micro-Fork Generator</a>.
</p>
<p>
Please see the <a href="https://github.com/hardcoded-fp/openpilot/">README for guidance and instructions</a>.
</p>
<div class="filter-bar">
<label for="filter-input"><strong>Filter vehicles</strong></label><br />
<input id="filter-input" type="search" placeholder="Type a model, make, or branch name" autocomplete="off" />
<div class="filter-meta"><span id="match-count"></span></div>
</div>
"""

    # Make it a nested list
    body = ""
    for base in base_cars_base_branches:
        escaped_base = html.escape(base)
        body += f'<section class="branch-section" data-branch="{escaped_base.lower()}">'
        body += f"<h2>{escaped_base}</h2>"
        body += "<ul>"
        sorted_cars = sorted(base_cars_base_branches[base].keys())
        for car in sorted_cars:
            escaped_car = html.escape(car)
            branch_name = base_cars_base_branches[base][car]
            escaped_branch_name = html.escape(branch_name)
            search_text = html.escape(f"{base} {car} {branch_name}".lower(), quote=True)
            body += f'<li class="car-item" data-search="{search_text}"><code>{escaped_car}</code>'
            body += f"<ul>"
            body += f"<li>Custom Software URL: <code>https://installer.comma.ai/hardcoded-fp/{escaped_branch_name}</code></li>"
            body += f'<li><a href="https://github.com/hardcoded-fp/openpilot/tree/{escaped_branch_name}">View on GitHub</a></li>'
            body += f"</ul>"
            body += f"</li>"

        body += "</ul>"
        body += "</section>"
    footer = """
<script>
const filterInput = document.getElementById("filter-input");
const matchCount = document.getElementById("match-count");
const sections = Array.from(document.querySelectorAll(".branch-section"));

function getQueryFromUrl() {
  const searchParams = new URLSearchParams(window.location.search);
  return searchParams.get("q") || "";
}

function updateSearchUrl(query) {
  const url = new URL(window.location.href);

  if (query === "") {
    url.searchParams.delete("q");
  } else {
    url.searchParams.set("q", query);
  }

  if (url.href !== window.location.href) {
    window.history.replaceState(null, "", url);
  }
}

function normalizeSearchText(value) {
  return value.trim().toLowerCase().replace(/[-_]+/g, " ").replace(/\\s+/g, " ");
}

function applyFilter(updateUrl = false) {
  const rawQuery = filterInput.value.trim();
  const query = normalizeSearchText(rawQuery);
  let totalVisible = 0;

  if (updateUrl) {
    updateSearchUrl(rawQuery);
  }

  for (const section of sections) {
    const items = Array.from(section.querySelectorAll(".car-item"));
    let visibleInSection = 0;

    for (const item of items) {
      const matches = query === "" || normalizeSearchText(item.dataset.search).includes(query);
      item.hidden = !matches;
      if (matches) {
        visibleInSection += 1;
      }
    }

    section.hidden = visibleInSection === 0;
    totalVisible += visibleInSection;
  }

  matchCount.textContent = query === ""
    ? `Showing all ${totalVisible} generated branches.`
    : `Showing ${totalVisible} matching generated branches.`;
}

filterInput.value = getQueryFromUrl();
filterInput.addEventListener("input", () => applyFilter(true));
window.addEventListener("popstate", () => {
  filterInput.value = getQueryFromUrl();
  applyFilter();
});
applyFilter();
</script>
</body>
</html>
"""
    # Make pages directory if it doesn't exist
    os.system("mkdir -p pages")
    with open("pages/index.html", "w") as f:
        f.write(header + body + footer)


def push_generated_branches(base_cars_base_branches):
    branch_names = []
    for branches_by_car in base_cars_base_branches.values():
        branch_names.extend(branches_by_car.values())

    unique_branch_names = sorted(set(branch_names))
    remote_tips = get_remote_branch_tips(unique_branch_names)
    local_tips = get_local_branch_tips(unique_branch_names)
    changed_branch_names = [
        branch_name
        for branch_name in unique_branch_names
        if local_tips[branch_name] != remote_tips.get(branch_name)
    ]
    skipped_count = len(unique_branch_names) - len(changed_branch_names)

    logging.info(
        "Skipping %s unchanged generated branches. Pushing %s changed branches.",
        skipped_count,
        len(changed_branch_names),
    )
    for batch in chunked(changed_branch_names, PUSH_BATCH_SIZE):
        logging.info("Pushing %s generated branches", len(batch))
        try:
            push_branch_batch_with_retries(batch)
        except subprocess.CalledProcessError:
            logging.warning(
                "Batch push failed. Falling back to per-branch pushes for this batch."
            )
            for branch_name in batch:
                logging.info("Pushing branch %s", branch_name)
                push_branch_with_retries(branch_name)


def main(push=True):
    resolved_branches = prepare_op_repo()
    base_branches = tuple(resolved_branches.keys())

    base_cars = {}
    for base in base_branches:
        base_cars[base] = parse_cars(resolved_branches[base])

    base_cars_base_branches = {}
    for base in base_branches:
        base_cars_base_branches[base] = {}
        for car in base_cars[base]:
            branch = generate_branch(base, resolved_branches[base], car)
            base_cars_base_branches[base][car] = branch
    logging.info("Done generating branches")

    # Log base_cars_base_branches
    logging.info("base_cars_base_branches:")
    logging.info(pprint.pformat(base_cars_base_branches))

    # Generate HTML output
    generate_html(base_cars_base_branches)

    if push:
        # Run the command to push to origin all the branches
        # Copy .git/config from this git repo to comma_openpilot repo
        # This might make GitHub Actions work
        run("cp .git/config comma_openpilot/.git/config")
        logging.info("Pushing branches to origin")
        push_generated_branches(base_cars_base_branches)


if __name__ == "__main__":
    # Check if args has dry run, if so, don't push
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--no-dry-run":
        main()
    else:
        main(push=False)
