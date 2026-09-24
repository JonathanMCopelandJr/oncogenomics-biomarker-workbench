"""Allow ``python -m onco_workbench`` as an alias for the ``obw`` command."""

from onco_workbench.cli import main

raise SystemExit(main())
