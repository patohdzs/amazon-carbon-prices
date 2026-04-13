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

echo "[preflight] Running macOS toolchain checks..."

check_cmd xcode-select "Run: xcode-select --install"
check_cmd xcrun "Run: xcode-select --install"
check_cmd clang "Install/reinstall Command Line Tools."
check_cmd make "Install Command Line Tools: xcode-select --install"
check_cmd brew "Install Homebrew from https://brew.sh/"
check_cmd Rscript "Install R: brew install r"
check_cmd gfortran "Install GNU Fortran: brew install gcc"

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

if [[ "$failures" -gt 0 ]]; then
    echo "[preflight] ${failures} check(s) failed."
    echo "[preflight] Fix the failures above, then rerun: ./scripts/preflight_macos.sh"
    exit 1
fi

echo "[preflight] All checks passed."
