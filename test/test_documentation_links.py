"""Guards against re-accumulating research-citation GitHub links in docs.

PlasmiCord's docs went through five rounds of "study these external repos" research,
each of which tended to leave behind markdown links to repos that are only cited or
accepted as opaque user-supplied evidence, not actually run or parsed by PlasmiCord's
own code. Those were all stripped back to plain text in one pass; this test keeps that
invariant from silently regressing in a future round -- see docs/REVIEW.md's five
"round of study" bullets for the removed context.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Only these repos are directly run by PlasmiCord's code or have their native output
# format parsed by it (mob-recon/tadrep native importers, the annotation.py engines),
# plus PlasmiCord's own repository. Every other GitHub link is a research citation and
# must stay plain text, not a markdown/clickable link.
ALLOWED_GITHUB_REPOS = {
    "phac-nml/mob-suite",
    "oschwengers/bakta",
    "tseemann/prokka",
    "ncbi/amr",
    "mikolmogorov/flye",
    "rrwick/unicycler",
    "oschwengers/tadrep",
    "ubeffiong/plasmicord",
}

SCANNED_SUFFIXES = {".md", ".py", ".js", ".toml", ".html", ".css", ".yml", ".yaml", ".json"}
EXCLUDED_DIR_NAMES = {".git", "build", "dist", "node_modules", "__pycache__"}
GITHUB_LINK = re.compile(r"github\.com/([A-Za-z0-9_.-]+)/([A-Za-z0-9_.-]+)")


def repo_files():
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in SCANNED_SUFFIXES:
            continue
        parts = set(path.relative_to(ROOT).parts)
        if parts & EXCLUDED_DIR_NAMES or any(p.startswith("results") or p.endswith(".egg-info") for p in parts):
            continue
        yield path


class DocumentationLinkTests(unittest.TestCase):
    def test_only_directly_integrated_tools_get_github_links(self):
        offenders = []
        for path in repo_files():
            text = path.read_text(encoding="utf-8", errors="ignore")
            for org, repo in GITHUB_LINK.findall(text):
                if repo.endswith(".git"):
                    repo = repo[:-len(".git")]
                slug = f"{org}/{repo}".lower()
                if slug not in ALLOWED_GITHUB_REPOS:
                    offenders.append(f"{path.relative_to(ROOT)}: github.com/{org}/{repo}")
        self.assertFalse(offenders,
            "Found GitHub links to repos PlasmiCord doesn't directly run/parse output from "
            "(research citations should be plain text, not links):\n" + "\n".join(offenders))


if __name__ == "__main__":
    unittest.main()
