## Commit Message Convention

We follow [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <description>

[optional body]

[optional footer(s)]
```

### Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, etc.)
- `refactor`: Code refactoring
- `perf`: Performance improvements
- `test`: Adding or updating tests
- `chore`: Maintenance tasks

### Examples:
```bash
git commit -m "feat: add support for DuckDB backend"
git commit -m "fix: resolve memory leak in scoring module"
git commit -m "docs: update installation instructions"
git commit -m "chore: update dependencies to latest versions"
```

### Breaking Changes:
```bash
git commit -m "feat!: remove deprecated API endpoints

BREAKING CHANGE: The /v1/score endpoint has been removed.
Use /v2/score instead."
```

## Dependencies and security checks

Use Python 3.11–3.13. Install the package and test tools in a virtual environment:

```bash
python -m pip install ".[all]" pytest pytest-regtest packaging
pytest -q
```

`pyproject.toml` declares supported dependencies and security minimums.
`requirements.txt` locks all extras for the Python 3.11 Linux release environment
(glibc 2.34 or newer); it is not a cross-platform or cross-Python lock file.
To refresh it using [uv](https://docs.astral.sh/uv/):

```bash
uv pip compile pyproject.toml --all-extras --python-version 3.11 \
  --python-platform x86_64-manylinux_2_34 --upgrade -o requirements.txt
python -m pip install pip-audit
pip-audit --strict --no-deps --disable-pip -r requirements.txt
```

Test dependency updates before committing the regenerated lock. CI tests
fresh resolutions on supported Python versions, audits the lock weekly, and
builds both distributions. Dependabot checks Python packages and pinned GitHub
Actions weekly. A clean dependency audit does not replace checking GitHub's
Dependabot, code-scanning, and secret-scanning alerts.

## Releasing to PyPI

Versions come from Git tags through `setuptools-scm`; do not maintain a separate
version in `pyproject.toml`. Source distributions retain their version when
rebuilt without Git. GitHub's automatically generated source archives are not
the supported source distributions; use the `.tar.gz` published on PyPI.

### One-time setup

The publishing workflow currently uses the existing repository secrets:

- `PYPI_USERNAME`: `__token__`
- `PYPI_PASSWORD`: a valid PyPI API token authorized to upload `easypqp`

PyPI account passwords are no longer supported. Keep tokens in GitHub Actions
secrets, never in source files or chat. Existing secrets can be reused only if
their token is still valid; account recovery may still be required if PyPI has
disabled the token or account. No Trusted Publisher or GitHub `pypi` environment
is required for this secret-based workflow. Attestations are disabled because
PyPI upload attestations require Trusted Publishing.

Docker publishing still requires `DOCKER_USERNAME` and `DOCKER_PASSWORD`
(use a scoped Docker Hub access token, not an account password).

### Optional later migration to Trusted Publishing

After PyPI account access is restored:

1. In the PyPI `easypqp` project's Publishing settings, add a **Trusted Publisher**:
   - Owner: `grosenberger`
   - Repository: `easypqp`
   - Workflow: `pythonpublish.yml`
   - Environment: `pypi`
2. Create the GitHub Actions environment `pypi`, restrict it to release tags,
   and configure required reviewers as appropriate.
3. Update the workflow's `publish` job to use that environment and grant
   `id-token: write`. Remove the credential-check step and the action's
   `user`, `password`, and `attestations: false` inputs to enable OIDC publishing
   and attestations.
4. After a successful Trusted Publishing upload, remove the obsolete
   `PYPI_USERNAME` / `PYPI_PASSWORD` secrets. Keep them until migration succeeds.

### Release procedure

1. Wait for CI to pass on the commit being released. Review `CHANGELOG.md`
   and include any desired updates in that commit; automation no longer
   pushes changelog commits directly onto `master`/`main`.
2. Push a new version tag, e.g. `0.1.62` or `v0.1.62`. Do not move existing tags.
3. The **Prepare draft release** workflow generates notes for that tag and
   creates a draft GitHub release. Review the notes and prerelease setting,
   then publish the draft as a maintainer. This human-triggered publication is
   intentional: publishing with `GITHUB_TOKEN` would not trigger downstream
   PyPI/Docker workflows.
4. **Upload Python Package** builds the exact tag on Python 3.11, checks wheel
   and sdist metadata against the tag, checks the packaged UniMod database,
   installs/tests the wheel with locked dependencies, and stores the validated
   distributions. A separate publishing job uploads those same files using
   the existing repository secrets. Docker runs independently.

For a failed upload, inspect the error before retrying. Manual dispatch of
**Upload Python Package** requires an explicit existing tag; it never silently
publishes the default branch. After migrating to Trusted Publishing, if the
`pypi` environment is restricted to tags, also select that tag in GitHub's
**Use workflow from** selector when dispatching.
Do not use `skip-existing` to hide version
mismatches. Existing tags predating the packaging fix cannot be repaired by
retagging; release a new version containing the fix instead.
