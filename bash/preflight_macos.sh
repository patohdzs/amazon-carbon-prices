#!/usr/bin/env bash

set -euo pipefail

if [[ "$(uname -s 2>/dev/null)" != "Darwin" ]]; then
    echo "[preflight] Non-macOS host detected; skipping macOS preflight."
    exit 0
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

failures=0

ok() {
    echo "[preflight] OK: $1"
}

warn() {
    echo "[preflight] WARN: $1"
}

fail() {
    echo "[preflight] FAIL: $1"
    failures=$((failures + 1))
}

check_cmd() {
    local cmd="$1"
    local hint="$2"
    if command -v "$cmd" >/dev/null 2>&1; then
        ok "$cmd found at $(command -v "$cmd")"
    else
        fail "$cmd not found. $hint"
    fi
}

required_r_version_from_project() {
    local version
    version=""

    if [[ -f ".R-version" ]]; then
        version="$(head -n 1 .R-version | tr -d '[:space:]')"
    elif [[ -f "renv.lock" ]]; then
        version="$(sed -n '1,40p' renv.lock \
            | sed -n 's/.*"Version"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' \
            | head -n 1)"
    fi

    if [[ -n "$version" ]]; then
        echo "$version"
        return 0
    fi
    return 1
}

echo "[preflight] Running macOS toolchain checks..."

check_cmd xcode-select "Run: xcode-select --install"
check_cmd xcrun "Run: xcode-select --install"
check_cmd clang "Install/reinstall Command Line Tools."
check_cmd make "Install Command Line Tools: xcode-select --install"
check_cmd brew "Install Homebrew from https://brew.sh/"
check_cmd rig "Install rig: brew install rig"
check_cmd pkg-config "Install pkg-config: brew install pkg-config"
check_cmd gfortran "Install GNU Fortran: brew install gcc"

required_r_version="$(required_r_version_from_project || true)"
if [[ -n "$required_r_version" ]]; then
    ok "Project requires R ${required_r_version} (from .R-version/renv.lock)."
else
    warn "Could not determine required R version from .R-version or renv.lock."
fi

if command -v Rscript >/dev/null 2>&1; then
    ok "Rscript found at $(command -v Rscript)"
    current_r_version="$(Rscript -e 'cat(paste0(R.version$major, ".", R.version$minor))' 2>/dev/null || true)"
    if [[ -n "$required_r_version" && "$current_r_version" != "$required_r_version" ]]; then
        fail "R version mismatch: current=${current_r_version}, required=${required_r_version}. Fix with: sudo rig add ${required_r_version} && sudo rig default ${required_r_version}"
    elif [[ -n "$current_r_version" ]]; then
        ok "R version is ${current_r_version}."
    else
        fail "Could not determine installed R version."
    fi
else
    if [[ -n "$required_r_version" ]]; then
        fail "Rscript not found. Install required R with: sudo rig add ${required_r_version} && sudo rig default ${required_r_version}"
    else
        fail "Rscript not found. Install R and rerun preflight."
    fi
fi

if [[ -n "${SDKROOT:-}" ]]; then
    fail "SDKROOT is set to '${SDKROOT}'. Unset it before bootstrap: unset SDKROOT"
fi
if [[ -n "${CPATH:-}" ]]; then
    fail "CPATH is set to '${CPATH}'. Unset it before bootstrap: unset CPATH"
fi
if [[ -n "${C_INCLUDE_PATH:-}" ]]; then
    fail "C_INCLUDE_PATH is set to '${C_INCLUDE_PATH}'. Unset it before bootstrap: unset C_INCLUDE_PATH"
fi
if [[ -n "${CPLUS_INCLUDE_PATH:-}" ]]; then
    fail "CPLUS_INCLUDE_PATH is set to '${CPLUS_INCLUDE_PATH}'. Unset it before bootstrap: unset CPLUS_INCLUDE_PATH"
fi
if [[ -n "${CPPFLAGS:-}" ]]; then
    warn "CPPFLAGS is set. If builds fail, unset it for bootstrap."
fi
if [[ -n "${CFLAGS:-}" ]]; then
    warn "CFLAGS is set. If builds fail, unset it for bootstrap."
fi
if [[ -n "${LDFLAGS:-}" ]]; then
    warn "LDFLAGS is set. If builds fail, unset it for bootstrap."
fi

if xcode_dir="$(xcode-select -p 2>/dev/null)"; then
    ok "xcode-select points to ${xcode_dir}"
else
    fail "xcode-select is not configured. Run: sudo xcode-select --switch /Library/Developer/CommandLineTools"
fi

sdk_path="$(xcrun --sdk macosx --show-sdk-path 2>/dev/null || true)"
if [[ -z "$sdk_path" ]]; then
    fail "Cannot resolve macOS SDK path with xcrun."
else
    ok "macOS SDK path: ${sdk_path}"
    if [[ -f "${sdk_path}/usr/include/string.h" ]]; then
        ok "C header check passed: ${sdk_path}/usr/include/string.h"
    else
        fail "Missing C headers in SDK (${sdk_path}/usr/include/string.h not found). Reinstall Command Line Tools."
    fi
fi

tmp_src="$(mktemp /tmp/preflight_c_XXXXXX.c)"
tmp_bin="$(mktemp /tmp/preflight_c_XXXXXX)"
cat > "$tmp_src" <<'C_SRC'
#include <string.h>
int main(void) { return (int)strlen("ok") == 2 ? 0 : 1; }
C_SRC

if clang "$tmp_src" -o "$tmp_bin" >/dev/null 2>&1; then
    ok "Minimal C compile test passed."
else
    fail "Minimal C compile test failed (headers/toolchain broken). Reinstall CLT."
fi
rm -f "$tmp_src" "$tmp_bin"

if [[ -f "Brewfile" ]]; then
    if brew bundle check --file Brewfile >/dev/null 2>&1; then
        ok "Brewfile dependencies are installed."
    else
        fail "Missing Homebrew dependencies from Brewfile. Run: brew bundle --file Brewfile"
    fi
else
    warn "No Brewfile found in repo root."
fi

im6_prefix="$(brew --prefix imagemagick@6 2>/dev/null || true)"
if [[ -z "$im6_prefix" ]]; then
    fail "imagemagick@6 not found. Install with: brew install imagemagick@6"
else
    im6_pc_dir="${im6_prefix}/lib/pkgconfig"
    im6_header="${im6_prefix}/include/ImageMagick-6/Magick++.h"
    if [[ ! -f "$im6_header" ]]; then
        fail "Missing Magick++ header at ${im6_header}. Reinstall imagemagick@6."
    else
        ok "Magick++ header found: ${im6_header}"
    fi
    if [[ -d "$im6_pc_dir" ]]; then
        PKG_CONFIG_PATH="${im6_pc_dir}:${PKG_CONFIG_PATH:-}" pkg-config --exists Magick++ \
            && ok "pkg-config can resolve Magick++." \
            || fail "pkg-config cannot resolve Magick++. Try: export PKG_CONFIG_PATH='${im6_pc_dir}:\$PKG_CONFIG_PATH'"
    else
        fail "Missing pkg-config directory at ${im6_pc_dir}."
    fi
fi

if [[ "$failures" -gt 0 ]]; then
    echo "[preflight] ${failures} check(s) failed."
    echo "[preflight] Fix the failures above, then rerun: ./bash/preflight_macos.sh"
    exit 1
fi

echo "[preflight] All checks passed."
