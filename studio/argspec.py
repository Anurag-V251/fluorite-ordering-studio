"""
argspec.py - describe the command-line options of a script as JSON, by running it up to its parse_args() call.

    python3 argspec.py SCRIPT.py [SUBCOMMAND]

The script's ArgumentParser is intercepted at parse_args / parse_known_args, so nothing of the script runs
beyond building its parser (module-level imports do run).  The web page builds its option forms from this, so
every script - present or future - gets a form without the app knowing its options.
"""
import argparse
import json
import os
import runpy
import sys


class _Captured(BaseException):      # BaseException: not swallowed by a script's `except Exception`
    def __init__(self, parser):
        self.parser = parser


def _intercept(self, *a, **k):
    raise _Captured(self)


def _jsonable(v):
    if v is None or isinstance(v, (bool, int, float, str)):
        return v
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    return str(v)


def _help(a, parser):
    h = a.help or ""
    if h is argparse.SUPPRESS:
        return ""
    if "%(" in h:
        try:
            return h % {**vars(a), "default": a.default, "prog": parser.prog}
        except (KeyError, ValueError, TypeError):
            return h
    return h.replace("%%", "%")


def describe_parser(parser):
    groups = {}
    for g in parser._action_groups:
        for a in g._group_actions:
            groups[id(a)] = g.title
    out = []
    for a in parser._actions:
        kind = type(a).__name__
        if kind == "_HelpAction" or kind == "_VersionAction":
            continue
        d = {
            "dest": a.dest,
            "flags": list(a.option_strings),
            "positional": not a.option_strings,
            "kind": {"_StoreAction": "store", "_StoreTrueAction": "true", "_StoreFalseAction": "false",
                     "_AppendAction": "append", "_CountAction": "count", "_StoreConstAction": "const",
                     "_SubParsersAction": "subcommand", "BooleanOptionalAction": "bool", "_ExtendAction": "extend",
                     "_AppendConstAction": "const"}.get(kind, kind),
            "nargs": a.nargs if not isinstance(a.nargs, int) else int(a.nargs),
            "type": getattr(a.type, "__name__", None) if a.type else None,
            "default": _jsonable(a.default) if a.default is not argparse.SUPPRESS else None,
            "choices": _jsonable(list(a.choices)) if a.choices is not None and kind != "_SubParsersAction" else None,
            "help": _help(a, parser),
            "required": bool(getattr(a, "required", False)) or (not a.option_strings and a.nargs not in ("?", "*")),
            "metavar": _jsonable(a.metavar),
            "group": groups.get(id(a)),
        }
        if kind == "_SubParsersAction":
            d["subcommands"] = {name: describe_parser(sp) for name, sp in a.choices.items()}
        out.append(d)
    return {"prog": parser.prog, "description": parser.description or "", "epilog": parser.epilog or "",
            "actions": out}


def main():
    script = os.path.abspath(sys.argv[1])
    sys.path.insert(0, os.path.dirname(script))
    argparse.ArgumentParser.parse_args = _intercept
    argparse.ArgumentParser.parse_known_args = _intercept
    sys.argv = [script]
    import io
    import contextlib
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            runpy.run_path(script, run_name="__main__")
    except _Captured as c:
        print(json.dumps({"ok": True, **describe_parser(c.parser)}))
        return
    except SystemExit as e:
        print(json.dumps({"ok": False, "error": f"script exited ({e.code}) before parsing arguments",
                          "output": buf.getvalue()[-2000:]}))
        return
    except Exception as e:     # noqa: BLE001
        print(json.dumps({"ok": False, "error": f"{type(e).__name__}: {e}", "output": buf.getvalue()[-2000:]}))
        return
    print(json.dumps({"ok": False, "error": "no argparse parser found", "output": buf.getvalue()[-2000:]}))


if __name__ == "__main__":
    main()
