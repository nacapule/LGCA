"""Put page-data.json into the page: island-study.html.

    python3 build/make_page.py

build/page.html is the page without its data. The result opens from a file, with no
server and no network.
"""
import json
import re
from pathlib import Path

STUDY = Path(__file__).resolve().parent.parent
SLOT = re.compile(r'(<script id="page-data" type="application/json">)(.*?)(</script>)', re.S)


def main():
    page = (STUDY / "build" / "page.html").read_text()
    data = (STUDY / "page-data.json").read_text()
    json.loads(data)
    # escape "<" so the data can never close the script tag; JSON reads it back the same
    data = data.replace("<", "\\u003c")
    page, n = SLOT.subn(lambda m: m.group(1) + data + m.group(3), page)
    if n != 1:
        raise SystemExit(f"build/page.html: expected one page-data slot, found {n}")
    target = STUDY / "island-study.html"
    target.write_text(page)
    print("wrote", target, f"{target.stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
