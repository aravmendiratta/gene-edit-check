# Releasing

How to publish a new version. Do the one-time setup once; after that a release is one tag.

## One-time setup

1. **Create the GitHub repository** and push this project to it.
2. **Replace the placeholder.** Search for `YOUR_GITHUB_USERNAME` and replace it with your GitHub username in
   `README.md`, `pyproject.toml`, `CITATION.cff` and `plugin.json`. The release workflow refuses to run
   while any placeholder is left.
3. **Make the figure visible.** The README shows `outputs/figure_crc-pd1-exhaustion.png` from the `main`
   branch, so it must be committed (it is) and pushed.
4. **Register on PyPI** at https://pypi.org, then go to *Publishing* and add a **pending trusted publisher**:
   - PyPI project name: `crispr-reflexion`
   - Owner: your GitHub username, repository: `crispr-reflexion`
   - Workflow name: `release.yml`
   - Environment name: `pypi`
5. In the GitHub repository, create an environment named `pypi` (Settings, then Environments).

Trusted publishing means no API token is stored anywhere.

## Each release

1. Update `version` in `pyproject.toml`, `plugin.json` and `skills/*/SKILL.md`, and `__version__` in
   `src/crispr_reflexion/__init__.py`. Add an entry to `CHANGELOG.md`.
2. Run `pytest`, then commit and push.
3. Tag and push the tag:

   ```bash
   git tag v0.3.0 && git push origin v0.3.0
   ```

4. The **Release** workflow builds the package, checks the tag matches the version, checks the reference
   table is inside the wheel, and publishes to PyPI.
5. Check it from a clean environment:

   ```bash
   pipx install crispr-reflexion && crispr-reflexion doctor
   ```

## Before the first public post

- Make the repository public and add a short description and topics (`alphagenome`, `crispr`,
  `single-cell`, `base-editing`).
- Open the live results page and share it from its Share menu if you want to link it.
- Do not quote the retracted figures (rho of +0.94 or +0.92, or "methylation-specific at PDCD1 and LAG3").
  See `docs/RESULTS.md` for what the data support.
