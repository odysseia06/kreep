#!/usr/bin/env bash
set -euo pipefail

if ! command -v gh >/dev/null 2>&1; then
  echo "gh CLI not found. Install from https://cli.github.com/"
  exit 1
fi

if ! gh auth status -h github.com >/dev/null 2>&1; then
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

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BACKLOG_FILE="${BACKLOG_FILE:-$ROOT_DIR/docs/triage/BACKLOG.md}"

if [[ ! -f "$BACKLOG_FILE" ]]; then
  echo "Backlog file not found: $BACKLOG_FILE"
  exit 1
fi

"$ROOT_DIR/scripts/triage_create_labels.sh" >/dev/null

mapfile -t titles < <(grep '^## ISSUE: ' "$BACKLOG_FILE" | sed 's/^## ISSUE: //')
if [[ ${#titles[@]} -eq 0 ]]; then
  echo "No issue sections found in $BACKLOG_FILE"
  exit 1
fi

existing_titles="$(gh issue list --repo "$REPO" --state all --limit 1000 --json title --jq '.[].title')"
existing_milestones="$(gh api "repos/$REPO/milestones?state=all&per_page=100" --paginate --jq '.[].title' 2>/dev/null || true)"

mapfile -t needed_milestones < <(grep '^Milestone: ' "$BACKLOG_FILE" | sed 's/^Milestone: //' | sort -u)
for m in "${needed_milestones[@]}"; do
  [[ -z "$m" ]] && continue
  if grep -Fxq "$m" <<<"$existing_milestones"; then
    continue
  fi
  gh api --method POST "repos/$REPO/milestones" -f title="$m" -f state="open" >/dev/null
  existing_milestones+=$'\n'"$m"
  echo "create milestone: $m"
done

extract_section() {
  local title="$1"
  awk -v t="$title" '
    $0 == "## ISSUE: " t { in_section=1; next }
    in_section && /^## ISSUE: / { exit }
    in_section { print }
  ' "$BACKLOG_FILE"
}

created=0
skipped=0

for title in "${titles[@]}"; do
  if grep -Fxq "$title" <<<"$existing_titles"; then
    echo "skip issue (exists): $title"
    skipped=$((skipped + 1))
    continue
  fi

  section="$(extract_section "$title")"
  if [[ -z "$section" ]]; then
    echo "skip issue (parse failed): $title"
    skipped=$((skipped + 1))
    continue
  fi

  labels_csv="$(printf '%s\n' "$section" | sed -n 's/^Labels: //p' | head -n 1)"
  milestone="$(printf '%s\n' "$section" | sed -n 's/^Milestone: //p' | head -n 1)"
  body="$(printf '%s\n' "$section" | sed -n '/^Evidence:/,$p')"

  if [[ -z "$body" ]]; then
    body="$section"
  fi

  tmp_body="$(mktemp)"
  {
    echo "Seeded from docs/triage/BACKLOG.md."
    echo
    printf '%s\n' "$body"
  } >"$tmp_body"

  args=(issue create --repo "$REPO" --title "$title" --body-file "$tmp_body")

  IFS=',' read -r -a label_parts <<<"$labels_csv"
  for raw in "${label_parts[@]}"; do
    label="$(sed 's/^ *//;s/ *$//' <<<"$raw")"
    [[ -z "$label" ]] && continue
    args+=(--label "$label")
  done

  if [[ -n "$milestone" ]]; then
    args+=(--milestone "$milestone")
  fi

  gh "${args[@]}" >/dev/null
  rm -f "$tmp_body"

  existing_titles+=$'\n'"$title"
  echo "create issue: $title"
  created=$((created + 1))
done

echo "issue sync complete for $REPO (created=$created, skipped=$skipped)"
