#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "--help" ]]; then
    printf 'Usage: bash %s [DESTINATION]\n' "$0"
    printf 'Recursively download all GLASS ET MODIS 1KM files, keeping YYYY/DDD folders.\n'
    printf 'Default destination: ./GLASS_ET_MODIS_1KM\n'
    printf 'Requires GNU Wget. Rerun with the same destination to resume.\n'
    exit 0
fi

if (( $# > 1 )); then
    printf 'Expected at most one destination argument. Use --help for usage.\n' >&2
    exit 2
fi

if ! command -v wget >/dev/null 2>&1; then
    printf 'GNU Wget is required. Install it in your Bash environment first.\n' >&2
    exit 1
fi

destination="${1:-./GLASS_ET_MODIS_1KM}"
mkdir -p -- "$destination"

printf 'Downloading all years and file types to: %s\n' "$destination"
printf 'Log: %s/download.log\n' "$destination"
printf 'Press Ctrl+C to stop; rerun the same command to resume.\n'

exec wget \
    --recursive \
    --level=0 \
    --no-parent \
    --no-host-directories \
    --cut-dirs=4 \
    --continue \
    --reject='index.html*' \
    --wait=1 \
    --random-wait \
    --tries=10 \
    --timeout=60 \
    --retry-connrefused \
    --directory-prefix="$destination" \
    --append-output="$destination/download.log" \
    'https://glass.hku.hk/archive/ET/MODIS/1KM/2008/033'