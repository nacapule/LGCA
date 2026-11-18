"""Build the standalone Island Study dashboard from the published data.

    python3 build/make_page.py

The script embeds page-data.json and entropy-page.json into build/page.html.
The data files must agree on the larger-box settings, seed counts and saved steps.
The result, island-study.html, opens directly in a browser.
"""
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
SLOT = re.compile(r'(<script id="page-data" type="application/json">)(.*?)(</script>)', re.S)
ENTROPY_SLOT = re.compile(r'(<script id="entropy-data" type="application/json">)(.*?)(</script>)', re.S)


def checkpoints(end):
    """Steps at which a run saved its lattice, through step `end`."""
    return [t for t in (0, 100, 500, 1000, 2000, 5000, 10000, 20000) if t <= end] + list(range(40000, end + 1, 20000))


def check_entropy(data, entropy):
    """Stop unless the entropy has every larger-box setting of the page data, with as many
    seeds and every saved step to the step its runs reached (so an entropy file older than
    the last merge is not used)."""
    have = {(e["alpha"], e["L"]): (len(e["seeds"]), e["steps"]) for e in entropy["entries"]}
    want = {(e["alpha"], e["L"]): (e["late"]["seeds"], checkpoints(e["horizon"])) for e in data["size"]["entries"]}
    if have != want:
        diff = sorted(k for k in set(have) | set(want) if have.get(k) != want.get(k))
        say = lambda v: f"{v[0]} seeds, {len(v[1])} saved steps to step {v[1][-1]}" if v else "none"
        raise SystemExit(f"entropy-page.json does not match page-data.json at {len(diff)} settings, "
                         f"such as alpha {diff[0][0]:g} L {diff[0][1]}: entropy {say(have.get(diff[0]))}, "
                         f"runs {say(want.get(diff[0]))}. Use matching data files.")


def main():
    page = (HERE / "build" / "page.html").read_text()
    data = (HERE / "page-data.json").read_text()
    entropy = (HERE / "entropy-page.json").read_text()
    check_entropy(json.loads(data), json.loads(entropy))
    # escape "<" so the data can never close the script tag; JSON reads it back the same
    data = data.replace("<", "\\u003c")
    page, n = SLOT.subn(lambda m: m.group(1) + data + m.group(3), page)
    if n != 1:
        raise SystemExit(f"build/page.html: expected one page-data slot, found {n}")
    entropy = entropy.replace("<", "\\u003c")
    page, n = ENTROPY_SLOT.subn(lambda m: m.group(1) + entropy + m.group(3), page)
    if n != 1:
        raise SystemExit(f"build/page.html: expected one entropy-data slot, found {n}")
    target = HERE / "island-study.html"
    target.write_text(page)
    print("wrote", target, f"{target.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
