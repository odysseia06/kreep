#!/usr/bin/env bash
set -euo pipefail

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=true
  echo "=== DRY RUN MODE (no changes will be made) ==="
fi

if ! command -v gh >/dev/null 2>&1; then
  echo "gh CLI not found. Install from https://cli.github.com/"
  exit 1
fi

if [[ "$DRY_RUN" == false ]] && ! gh auth status -h github.com >/dev/null 2>&1; then
  echo "gh is not authenticated for github.com. Run: gh auth login"
  exit 1
fi

REPO="${GH_REPO:-}"
if [[ -z "$REPO" ]]; then
  REPO="$(gh repo view --json nameWithOwner --jq '.nameWithOwner')"
fi

if [[ -z "$REPO" ]]; then
  echo "Could not resolve repository. Set GH_REPO=owner/name and retry."
  exit 1
fi

existing_labels=""
if [[ "$DRY_RUN" == false ]]; then
  existing_labels="$(gh label list --repo "$REPO" --limit 1000 --json name --jq '.[].name')"
fi

create_if_missing() {
  local name="$1"
  local color="$2"
  local desc="$3"

  if [[ -n "$existing_labels" ]] && grep -Fxq "$name" <<<"$existing_labels"; then
    echo "skip label: $name"
    return 0
  fi

  if [[ "$DRY_RUN" == true ]]; then
    echo "would create label: $name ($color) — $desc"
    return 0
  fi

  gh label create "$name" \
    --repo "$REPO" \
    --color "$color" \
    --description "$desc" >/dev/null

  existing_labels+=$'\n'"$name"
  echo "create label: $name"
}

while IFS='|' read -r name color desc; do
  [[ -z "$name" ]] && continue
  create_if_missing "$name" "$color" "$desc"
done <<'LABELS'
bug|d73a4a|Incorrect behavior, failing tests, or regressions
feature|0e8a16|New functionality or capability
perf|fbca04|Performance optimization or regression work
tech-debt|6f42c1|Maintainability and engineering quality work
docs|0075ca|Documentation or example improvements
security|b60205|Side-channel, misuse-hazard, or cryptographic safety work
P0-correctness/blocker|b60205|Must-fix correctness issue or workflow blocker
P1-important|fbca04|Important near-term work
P2-later|0e8a16|Lower-priority or roadmap work
area:traits-algebra|5319e7|Algebraic traits and abstractions
area:field-fp|1d76db|Prime field implementation (Fp)
area:ext-fields|0e8a16|Extension fields and tower fields
area:polynomials|aa6708|Polynomial algorithms and APIs
area:ntt-fft|0052cc|NTT/FFT-related functionality
area:factorization|c2e0c6|Polynomial factorization and roots
area:discrete-log|d4c5f9|Discrete logarithm functionality
area:serialization|006b75|Serde and data representation boundaries
area:randomness|f9d0c4|Random sampling and RNG usage
area:constant-time|b60205|Constant-time behavior and side-channel hardening
area:no-std|0366d6|no_std/alloc portability and feature gating
area:benchmarks|1d76db|Criterion/benchmark coverage and performance telemetry
area:docs-examples|0075ca|README/rustdoc/examples onboarding quality
area:infra-ci|6f42c1|Tooling, CI, workflows, and release automation
status:needs-investigation|e99695|Needs confirmation, reproducer, or scope refinement
status:blocked|000000|Blocked by dependency or prerequisite decision
LABELS

echo "label sync complete for $REPO"
