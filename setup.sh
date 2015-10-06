# This should be sourced after setting up a base release.

# Add the test release to PATH, LD_LIBRARY_PATH and PYTHON_PATH
source $(dirname $BASH_SOURCE)/bin/setup_mu2e_project.sh

# Define the test release environment variables.
source  ${MU2E_BASE_RELEASE}/bin/addlocal.sh
