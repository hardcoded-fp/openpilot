#!/usr/bin/env python3

import ast
import datetime
import logging
import os
import pprint

logging.basicConfig(level=logging.INFO, format="%(message)s")


def parse_cars(branch):
    """
    Parse the file with ast for a class named "CAR" and get a list of all the cars.
    We are looking for the values in quotes.

    Exerpt:

    class CAR:
      # Hyundai
      ELANTRA = "HYUNDAI ELANTRA 2017"
      ELANTRA_2021 = "HYUNDAI ELANTRA 2021"
      ELANTRA_HEV_2021 = "HYUNDAI ELANTRA HYBRID 2021"
      HYUNDAI_GENESIS = "HYUNDAI GENESIS 2015-2016"
      IONIQ = "HYUNDAI IONIQ HYBRID 2017-2019"

    Another format is to look in in fingerprints.py instead of values.py

    Exerpt:

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
    os.system(f"cd comma_openpilot && git checkout --force {branch}")

    # Get a list of values.py underneath the folder
    # "comma_openpilot/selfdrive/car/"

    values_py_paths = []
    fingerprints_py_paths = []
    cars = []

    for root, dirs, files in os.walk("comma_openpilot/selfdrive/car/"):
        values_py_paths += [os.path.join(root, f) for f in files if f == "values.py"]


    for path in values_py_paths:
        logging.info("Parsing %s", path)
        with open(path, "r") as f:
            tree = ast.parse(f.read())
            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef) and node.name == "CAR":
                    for c in node.body:
                        if isinstance(c, ast.Assign):
                            if isinstance(c.value, ast.Str):
                                cars.append(c.value.s)
                            # Sometimes it's an object initializer,
                            # If so, use the first argument
                            elif isinstance(c.value, ast.Call):
                                # Sometimes
                                if len(c.value.args) > 0 and isinstance(c.value.args[0], ast.Str):
                                    cars.append(c.value.args[0].s)

    for root, dirs, files in os.walk("comma_openpilot/selfdrive/car/"):
        fingerprints_py_paths += [os.path.join(root, f) for f in files if f == "fingerprints.py"]

    for path in fingerprints_py_paths:
        logging.info("Parsing %s", path)
        with open(path, "r") as f:
            tree = ast.parse(f.read())
            for node in ast.walk(tree):
                if isinstance(node, ast.Assign) and isinstance(node.value, ast.Dict):
                    for key in node.value.keys:
                        if isinstance(key, ast.Attribute):
                            cars.append(key.attr)

    # Log the cars
    logging.info("Found %d cars in %s", len(cars), branch)

    return cars


def prepare_op_repo():
    """
    Prepare the openpilot repo with master-ci and release3 branches
    """
    # Try to clone the repo to comma_openpilot.
    # If it fails, it means it already exists, so we can just pull
    # the latest changes.
    logging.info("Setting up openpilot repo. Ignore errors if it already exists.")

    os.system(
        "git clone -b master-ci https://github.com/commaai/openpilot.git comma_openpilot"
    )
    # Make sure that comma_openpilot is usiing that as the origin.
    os.system(
        "cd comma_openpilot && git remote set-url origin https://github.com/commaai/openpilot.git"
    )
    os.system("cd comma_openpilot && git fetch origin")
    os.system(
        "cd comma_openpilot && git checkout release3 && git reset --hard origin/release3"
    )
    os.system(
        "cd comma_openpilot && git checkout master-ci && git reset --hard origin/master-ci"
    )

    logging.info("Done setting up openpilot repo.")


def generate_branch(base, car):
    """
    Make a new branch for the car with a hardcoded fingerprint
    """

    # - instead of _ because the keyboard is one tap for - vs two for _
    # & is AND because & may be too special
    # Lowercase because there's no caps lock in the keyboard
    # Remove () because they are special characters and may cause issues
    branch_name = f"{base}-{car.replace(' ', '-').replace('&', 'AND').replace('(', '').replace(')','').replace('_', '-').lower()}"
    logging.info("Generating branch %s", branch_name)
    # Delete branch if it already exists
    os.system(f"cd comma_openpilot && git branch -D {branch_name}")
    # Make branch off of base branch
    os.system(
        f"cd comma_openpilot && git checkout {base} && git checkout -b {branch_name}"
    )
    # Make sure base branch is clean
    os.system(f"cd comma_openpilot && git reset --hard origin/{base}")
    # Append 'export FINGERPRINT="car name"' to the end of launch_env.sh
    os.system(f"echo 'export FINGERPRINT=\"{car}\"' >> comma_openpilot/launch_env.sh")
    # Commit the changes
    # Get date of current commit
    commit_date = os.popen(
        "cd comma_openpilot && git log -1 --format=%cd --date=iso-strict"
    ).read()
    author_date = os.popen(
        "cd comma_openpilot && git log -1 --format=%ad --date=iso-strict"
    ).read()

    os.system(
        f"cd comma_openpilot && git add launch_env.sh && GIT_AUTHOR_DATE='{author_date}' GIT_COMMITTER_DATE='{commit_date}' git commit -m 'Hardcode fingerprint for {car}'"
    )
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
"""

    # Make it a nested list
    body = ""
    for base in base_cars_base_branches:
        body += f"<h2>{base}</h2>"
        body += "<ul>"
        sorted_cars = sorted(base_cars_base_branches[base].keys())
        for car in sorted_cars:
            body += f"<li><code>{car}</code>"
            body += f"<ul>"
            body += f"<li>Custom Software URL: <code>https://installer.comma.ai/hardcoded-fp/{base_cars_base_branches[base][car]}</code></li>"
            body += f'<li><a href="https://github.com/hardcoded-fp/openpilot/tree/{base_cars_base_branches[base][car]}">View on GitHub</a></li>'
            body += f"</ul>"
            body += f"</li>"

        body += "</ul>"
    footer = """
</body>
</html>
"""
    # Make pages directory if it doesn't exist
    os.system("mkdir -p pages")
    with open("pages/index.html", "w") as f:
        f.write(header + body + footer)


def main(push=True):
    prepare_op_repo()

    base_cars = {}
    base_branches = ["master-ci", "nightly", "release3"]
    for base in base_branches:
        base_cars[base] = parse_cars(base)

    base_cars_base_branches = {}
    for base in base_branches:
        base_cars_base_branches[base] = {}
        for car in base_cars[base]:
            branch = generate_branch(base, car)
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
        os.system("cp .git/config comma_openpilot/.git/config")
        logging.info("Pushing branches to origin")
        os.system("cd comma_openpilot && git push origin --force --all")


if __name__ == "__main__":
    # Check if args has dry run, if so, don't push
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--no-dry-run":
        main()
    else:
        main(push=False)
