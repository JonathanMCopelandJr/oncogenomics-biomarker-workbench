# Security policy

## Scope

This is a local research and education tool that works on **synthetic** data:

- It needs **no credentials**, API keys, or network access. It downloads nothing at
  runtime.
- The Streamlit dashboard binds to `localhost` only, has no authentication or upload
  feature, and turns off Streamlit usage statistics.
- It is **not** designed or tested for deployment on a server, for multi-user use, or for
  processing real patient data or protected health information. Please do not use it
  that way.

## Supported versions

Only the latest commit on the `main` branch is supported. No versions have been released yet.

## Reporting a vulnerability

Please **do not open a public issue** for security problems. Instead:

1. Use GitHub's **private vulnerability reporting** ("Security" tab, then "Report a
   vulnerability") if it is enabled for this repository, or
2. contact the maintainer privately through their GitHub profile.

Include a description, steps to reproduce, and the affected version or commit. **Never
include real patient data, credentials, or other sensitive information in a report.**
Reports are acknowledged as soon as possible. Fixes are prioritized by severity.

## Data safety reminders

- Keep real or restricted data out of this repository. `data/raw/`, `data/processed/`,
  and `outputs/` are ignored by Git for this reason.
- Never commit `.env` files, `.streamlit/secrets.toml`, keys, or tokens. They are also
  ignored by Git.
