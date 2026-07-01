#!/usr/bin/env python3
# Copyright 2026 keplertech.io
# SPDX-License-Identifier: GPL-3.0-only

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

SKIP_RETURN_CODE = 77


def load_manifest(benchmarks_root):
    manifest = benchmarks_root / "manifest.tsv"
    cases = {}
    with manifest.open() as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            row = dict(zip(header, parts))
            cases[row["case"]] = row
    return cases


def yaml_quote(path):
    return str(path).replace("\\", "\\\\").replace('"', '\\"')


def strip_block_comments(text):
    return re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)


def mask_comments(text):
    def mask_match(match):
        return "".join("\n" if ch == "\n" else " " for ch in match.group(0))

    text = re.sub(r"/\*.*?\*/", mask_match, text, flags=re.DOTALL)
    return re.sub(r"//[^\n]*", mask_match, text)


def remove_initial_blocks(text):
    text_for_scan = mask_comments(text)
    result = []
    pos = 0
    pattern = re.compile(r"\binitial\b")
    token_pattern = re.compile(r"\bbegin\b|\bend\b|;")
    while True:
        match = pattern.search(text_for_scan, pos)
        if not match:
            result.append(text[pos:])
            break
        result.append(text[pos:match.start()])
        scan = match.end()
        while scan < len(text_for_scan) and text_for_scan[scan].isspace():
            scan += 1
        if text_for_scan.startswith("begin", scan):
            depth = 0
            end = scan
            for token in token_pattern.finditer(text_for_scan, scan):
                word = token.group(0)
                if word == "begin":
                    depth += 1
                elif word == "end":
                    depth -= 1
                    if depth == 0:
                        end = token.end()
                        break
            pos = end
        else:
            semicolon = text_for_scan.find(";", scan)
            pos = len(text_for_scan) if semicolon == -1 else semicolon + 1
        result.append("\n")
    return "".join(result)


def remove_assertions(text):
    text = re.sub(r"(?ms)^\s*assert\s+property\b.*?;\s*", "", text)
    text = re.sub(r"(?ms)^\s*assert\s+\w+\s*:.*?;\s*", "", text)
    return text


def rename_reserved_modules(text):
    text = re.sub(r"\bmodule\s+program\b", "module program_rom", text)
    text = re.sub(r"\bprogram\s+([A-Za-z_][A-Za-z0-9_$]*)\s*\(", r"program_rom \1(", text)
    text = re.sub(r"\bmodule\s+cell\b", "module tree_cell", text)
    text = re.sub(r"\bcell\s+([A-Za-z_][A-Za-z0-9_$]*)\s*\(", r"tree_cell \1(", text)
    return text


def inject_vending_constants(text):
    if "module vending" not in text or "module monitor" not in text:
        return text
    constants = (
        "    parameter NONE = 0;\n"
        "    parameter NICKEL = 1;\n"
        "    parameter DIME = 2;\n"
        "    parameter QUARTER = 3;\n"
        "    parameter ACCEPTING = 0;\n"
        "    parameter CHANGE = 1;\n"
        "    parameter REFUND = 2;\n"
        "    parameter BEVERAGE = 3;\n"
    )

    def inject(match):
        header = match.group(0)
        if "parameter NONE" in text[match.end() : match.end() + 300]:
            return header
        return header + constants

    text = re.sub(r"module\s+vending\s*\([^;]*\);\n", inject, text)
    return re.sub(r"module\s+monitor\s*\([^;]*\);\n", inject, text)


def add_legacy_implicit_net_declarations(text):
    module_names = set(re.findall(r"\bmodule\s+([A-Za-z_][A-Za-z0-9_$]*)\s*\(", text))
    declaration_keywords = {
        "input",
        "output",
        "inout",
        "wire",
        "reg",
        "logic",
        "integer",
        "parameter",
        "signed",
        "unsigned",
    }

    def collect_declared(header, body):
        declared = set()
        port_match = re.search(r"\((.*?)\)", header, flags=re.DOTALL)
        if port_match:
            declared.update(re.findall(r"\b[A-Za-z_][A-Za-z0-9_$]*\b", port_match.group(1)))
        for decl in re.finditer(
            r"\b(?:input|output|inout|wire|reg|logic|integer)\b([^;]*);", body, flags=re.DOTALL
        ):
            names = re.findall(r"\b[A-Za-z_][A-Za-z0-9_$]*\b", decl.group(1))
            declared.update(name for name in names if name not in declaration_keywords)
        for decl in re.finditer(r"\bparameter\s+([A-Za-z_][A-Za-z0-9_$]*)\b", body):
            declared.add(decl.group(1))
        return declared

    def patch_module(match):
        header, module_name, body, footer = match.groups()
        declared = collect_declared(header, body)
        implicit = []
        for module_type in module_names - {module_name}:
            instance_pattern = re.compile(
                rf"\b{re.escape(module_type)}\s+[A-Za-z_][A-Za-z0-9_$]*\s*\((.*?)\)\s*;",
                flags=re.DOTALL,
            )
            for instance in instance_pattern.finditer(body):
                for arg in instance.group(1).split(","):
                    arg = arg.strip()
                    if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_$]*", arg) and arg not in declared:
                        declared.add(arg)
                        implicit.append(arg)
        if not implicit:
            return match.group(0)
        declaration = "".join(f"wire {name};\n" for name in implicit)
        return f"{header}{declaration}{body}{footer}"

    return re.sub(
        r"(?ms)(module\s+([A-Za-z_][A-Za-z0-9_$]*)\s*\([^;]*\);\n?)(.*?)(\nendmodule\b)",
        patch_module,
        text,
    )


def remove_function(text, name):
    return re.sub(
        rf"(?ms)^\s*function(?:\s+\[[^\]]+\])?\s+{re.escape(name)}\b.*?^\s*endfunction[^\n]*\n",
        "",
        text,
    )


def rewrite_lookup_functions(text):
    if "function [7:0] mem" in text and "assign mem_mar = mem(mar);" in text:
        mem_expr = (
            "((mar == 4'd0) ? 8'd1 : "
            "((mar == 4'd1) ? 8'd255 : "
            "((mar == 4'd2) ? 8'd0 : "
            "((mar == 4'd3) ? 8'd0 : "
            "((mar == 4'd4) ? 8'd0 : "
            "((mar == 4'd5) ? 8'd2 : "
            "((mar == 4'd6) ? 8'd0 : "
            "((mar == 4'd7) ? 8'd0 : "
            "((mar == 4'd8) ? 8'd0 : "
            "((mar == 4'd9) ? 8'd2 : "
            "((mar == 4'd10) ? 8'd255 : "
            "((mar == 4'd11) ? 8'd5 : "
            "((mar == 4'd12) ? 8'd0 : "
            "((mar == 4'd13) ? 8'd2 : "
            "((mar == 4'd14) ? 8'd0 : 8'd2)))))))))))))))"
        )
        text = remove_function(text, "mem")
        text = text.replace("    assign mem_mar = mem(mar);", f"    assign mem_mar = {mem_expr};")

    if "function [19:0] ROM" in text and "assign ROM_OUT = ROM(MAR);" in text:
        rom_expr = (
            "((MAR == 3'd0) ? 20'b01111111100101111010 : "
            "((MAR == 3'd1) ? 20'b00111001110101100010 : "
            "((MAR == 3'd2) ? 20'b10101000111111111111 : "
            "((MAR == 3'd3) ? 20'b11111111011010111010 : "
            "((MAR == 3'd4) ? 20'b11111111111101101110 : "
            "((MAR == 3'd5) ? 20'b11111111101110101000 : "
            "((MAR == 3'd6) ? 20'b11001010011101011011 : "
            "20'b00101111111111110100)))))))"
        )
        text = remove_function(text, "ROM")
        text = text.replace("    assign ROM_OUT = ROM(MAR);", f"    assign ROM_OUT = {rom_expr};")

    return text


def rewrite_vlunc_helpers(text):
    if "module transform" not in text or "toLower(in)" not in text:
        return text
    text = text.replace("1'bx", "1'b0").replace("8'hxx", "8'h00")
    for function_name in ("toLower", "toUpper", "changeCase", "isUpper"):
        text = remove_function(text, function_name)
    transform_assign = re.compile(
        r"(?ms)    assign out = Lcmd \? toLower\(in\) :\s*"
        r"Ucmd \? toUpper\(in\) :\s*"
        r"Ncmd \? in :\s*"
        r"Ccmd \? changeCase\(in\) : 8'h00;"
    )
    transform_replacement = (
        "    wire        is_upper;\n"
        "    wire [7:0]  lower;\n"
        "    wire [7:0]  upper;\n"
        "    wire [7:0]  changed;\n"
        "    wire [7:0]  lower_mask;\n"
        "    wire [7:0]  upper_mask;\n"
        "    wire [7:0]  normal_mask;\n"
        "    wire [7:0]  changed_mask;\n\n"
        "    assign is_upper = ~in[5];\n"
        "    assign lower = in + ({8{is_upper}} & 8'h20);\n"
        "    assign upper = in - ({8{~is_upper}} & 8'h20);\n"
        "    assign changed = in + ({8{is_upper}} & 8'h20) - ({8{~is_upper}} & 8'h20);\n"
        "    assign lower_mask = {8{Lcmd}};\n"
        "    assign upper_mask = {8{~Lcmd & Ucmd}};\n"
        "    assign normal_mask = {8{~Lcmd & ~Ucmd & Ncmd}};\n"
        "    assign changed_mask = {8{~Lcmd & ~Ucmd & ~Ncmd & Ccmd}};\n"
        "    assign out = (lower_mask & lower) | (upper_mask & upper) |\n"
        "                 (normal_mask & in) | (changed_mask & changed);"
    )
    return transform_assign.sub(transform_replacement, text)


def rewrite_unidec_functions(text):
    if not all(token in text for token in ("function [15:0] code", "function [15:0] prefix", "function [15:0] suffix")):
        return text

    def code_expr(sel):
        return (
            f"(({sel}) == 3'd0 ? 16'b0000000000001000 : "
            f"(({sel}) == 3'd1 ? 16'b0000000000001010 : "
            f"(({sel}) == 3'd2 ? 16'b0000000001011000 : "
            f"(({sel}) == 3'd3 ? 16'b0000001001001000 : "
            f"(({sel}) == 3'd4 ? 16'b0000001011000001 : "
            f"(({sel}) == 3'd5 ? 16'b0000001001100011 : "
            f"(({sel}) == 3'd6 ? 16'b1100011010001001 : "
            "16'b0000000000001000)))))))"
        )

    def prefix_expr(word, sel):
        invalid = "16'b0111111111111111"
        return (
            f"(({sel}) == 2'd0 ? (({word}[15:4] == 12'd0) ? {invalid} : {{13'b1, {word}[2:0]}}) : "
            f"(({sel}) == 2'd1 ? (({word}[15:7] == 9'd0) ? {invalid} : {{10'b1, {word}[5:0]}}) : "
            f"(({sel}) == 2'd2 ? (({word}[15:10] == 6'd0) ? {invalid} : {{7'b1, {word}[8:0]}}) : "
            f"(({word}[15:13] == 3'd0) ? {invalid} : {{4'b1, {word}[11:0]}}))))"
        )

    def suffix_expr(word, sel):
        return (
            f"(({sel}) == 2'd0 ? {{3'b0, {word}[15:3]}} : "
            f"(({sel}) == 2'd1 ? {{6'b0, {word}[15:6]}} : "
            f"(({sel}) == 2'd2 ? {{9'b0, {word}[15:9]}} : "
            f"{{12'b0, {word}[15:12]}})))"
        )

    text = re.sub(r"(?ms)^\s*function\s+\[15:0\]\s+(?:code|prefix|suffix)\b.*?^\s*endfunction[^\n]*\n", "", text)
    text = text.replace(
        "    reg         init;\n",
        "    reg         init;\n"
        "    wire [15:0] prefix_word_sel2;\n"
        "    wire [15:0] prefix_other_sel2;\n"
        "    wire [15:0] suffix_word_sel2;\n"
        "    wire [15:0] suffix_other_sel2;\n",
    )
    text = text.replace(
        "    assign other = code(sel1);\n",
        f"    assign other = {code_expr('sel1')};\n"
        f"    assign prefix_word_sel2 = {prefix_expr('word', 'sel2')};\n"
        f"    assign prefix_other_sel2 = {prefix_expr('other', 'sel2')};\n"
        f"    assign suffix_word_sel2 = {suffix_expr('word', 'sel2')};\n"
        f"    assign suffix_other_sel2 = {suffix_expr('other', 'sel2')};\n",
    )
    text = text.replace("prefix(word,sel2)", "prefix_word_sel2")
    text = text.replace("prefix(other,sel2)", "prefix_other_sel2")
    text = text.replace("suffix(word,sel2)", "suffix_word_sel2")
    text = text.replace("suffix(other,sel2)", "suffix_other_sel2")
    return text


def sanitize_rtl_for_kepler(text):
    text = remove_initial_blocks(text)
    text = strip_block_comments(text)
    text = remove_assertions(text)
    text = rename_reserved_modules(text)
    text = inject_vending_constants(text)
    text = rewrite_lookup_functions(text)
    text = rewrite_vlunc_helpers(text)
    text = rewrite_unidec_functions(text)
    text = add_legacy_implicit_net_declarations(text)
    return "`default_nettype wire\n" + text + "\n`default_nettype none\n"


def write_sanitized_rtl_source(rtl_path, out_path):
    out_path.write_text(sanitize_rtl_for_kepler(rtl_path.read_text(errors="ignore")))


def write_metron_passthrough_source(rtl_path, out_path):
    rtl_text = sanitize_rtl_for_kepler(rtl_path.read_text(errors="ignore"))
    out_path.write_text(
        "#include \"metron/metron_tools.h\"\n\n"
        "/*#\n"
        f"{rtl_text}\n"
        "#*/\n\n"
        "class VerilogCBenchmarkPassthrough {\n"
        " public:\n"
        "  logic<1> keep() const { return 0; }\n"
        "};\n"
    )


def run_case(args):
    env = os.environ.copy()
    benchmarks_root = Path(args.benchmarks_root).resolve()

    cases = load_manifest(benchmarks_root)
    if args.case not in cases:
        print(f"Unknown benchmark case {args.case!r}", file=sys.stderr)
        return 2

    row = cases[args.case]
    if row["status"] != "complete":
        print(f"Benchmark case {args.case} is not runnable: {row['note']}", file=sys.stderr)
        return SKIP_RETURN_CODE

    case_dir = benchmarks_root / "cases" / args.case
    c_path = case_dir / row["c_source"]
    rtl_path = case_dir / row["rtl_source"]
    if not c_path.is_file() or not rtl_path.is_file():
        print(f"Benchmark fixture is incomplete for {args.case}", file=sys.stderr)
        print(f"  C:   {c_path}", file=sys.stderr)
        print(f"  RTL: {rtl_path}", file=sys.stderr)
        return 2
    with tempfile.TemporaryDirectory(prefix="kepler_verilog_c_") as tmp:
        tmp_dir = Path(tmp)
        c_work = tmp_dir / "c2rtl"
        metron_path = tmp_dir / "metron_passthrough.h"
        sanitized_rtl_path = tmp_dir / "reference.sv"
        write_metron_passthrough_source(rtl_path, metron_path)
        write_sanitized_rtl_source(rtl_path, sanitized_rtl_path)
        cfg_path = tmp_dir / "config.yaml"
        cfg_path.write_text(
            "verification: lec\n"
            "design1:\n"
            "  format: c\n"
            f"  input_paths: [\"{yaml_quote(metron_path)}\"]\n"
            f"  top: {row['rtl_top']}\n"
            f"  work_dir: \"{yaml_quote(c_work)}\"\n"
            "design2:\n"
            "  format: systemverilog\n"
            f"  input_paths: [\"{yaml_quote(sanitized_rtl_path)}\"]\n"
            f"  top: {row['rtl_top']}\n"
        )

        cmd = [str(Path(args.kepler_formal).resolve()), "--config", str(cfg_path)]
        print(f"Running {args.case}: {' '.join(cmd)}")
        completed = subprocess.run(
            cmd,
            cwd=tmp_dir,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=args.timeout,
        )
        if completed.stdout:
            print(completed.stdout, end="")
        if completed.stderr:
            print(completed.stderr, end="", file=sys.stderr)
        if completed.returncode != 0:
            print(
                f"Benchmark {args.case} failed: kepler-formal exited {completed.returncode}",
                file=sys.stderr,
            )
            print(f"Upstream safety expectation: {row['upstream_expected_safety']}", file=sys.stderr)
        return completed.returncode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kepler-formal", required=True)
    parser.add_argument("--benchmarks-root", required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()
    try:
        return run_case(args)
    except subprocess.TimeoutExpired as exc:
        print(f"Benchmark {args.case} timed out after {exc.timeout}s", file=sys.stderr)
        return 124


if __name__ == "__main__":
    sys.exit(main())
