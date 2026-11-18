#!/usr/bin/env bash
# Install the Quarto CLI into ~/.local/share/quarto and link ~/.local/bin/quarto.
# Kept outside ~/research so OneDrive never touches binary toolchains or symlinks.

set -euo pipefail

QUARTO_VERSION="${QUARTO_VERSION:-1.11.9}"
INSTALL_DIR="${HOME}/.local/share/quarto-${QUARTO_VERSION}"
BIN_DIR="${HOME}/.local/bin"

if command -v quarto >/dev/null 2>&1; then
  echo "Quarto already on PATH: $(command -v quarto) ($(quarto --version))"
  exit 0
fi

if [[ -x "${BIN_DIR}/quarto" ]]; then
  echo "Quarto already installed at ${BIN_DIR}/quarto ($("${BIN_DIR}/quarto" --version))"
  exit 0
fi

mkdir -p "${HOME}/.local/share" "${BIN_DIR}"

OS="$(uname -s)"
if [[ "${OS}" == "Darwin" ]]; then
  TARBALL="quarto-${QUARTO_VERSION}-macos.tar.gz"
elif [[ "${OS}" == "Linux" ]]; then
  ARCH="$(uname -m)"
  if [[ "${ARCH}" == "aarch64" || "${ARCH}" == "arm64" ]]; then
    TARBALL="quarto-${QUARTO_VERSION}-linux-arm64.tar.gz"
  else
    TARBALL="quarto-${QUARTO_VERSION}-linux-amd64.tar.gz"
  fi
else
  echo "Unsupported OS: ${OS}" >&2
  exit 1
fi

URL="https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/${TARBALL}"
TMP_TAR="$(mktemp)"
trap 'rm -f "${TMP_TAR}"' EXIT

echo "Downloading Quarto ${QUARTO_VERSION} from ${URL} ..."
curl -fsSL "${URL}" -o "${TMP_TAR}"

rm -rf "${INSTALL_DIR}"
mkdir -p "${INSTALL_DIR}"
tar -xzf "${TMP_TAR}" -C "${INSTALL_DIR}" --strip-components=1

# Remove macOS quarantine attribute if present so non-interactive builds work cleanly
if [[ "${OS}" == "Darwin" ]] && command -v xattr >/dev/null 2>&1; then
  xattr -dr com.apple.quarantine "${INSTALL_DIR}" 2>/dev/null || true
fi

ln -sf "${INSTALL_DIR}/bin/quarto" "${BIN_DIR}/quarto"
echo "Installed Quarto $("${BIN_DIR}/quarto" --version) -> ${BIN_DIR}/quarto"
