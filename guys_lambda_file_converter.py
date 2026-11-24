#!/usr/bin/env python3
"""Convert AMBER pmemd CpHMD lambda trajectories: .lambda <-> .lmb (LMB1).

LMB1 has eight newline-terminated ASCII header lines, then uint8 lambda
frames. Each frame uses ntitr bytes; steps are step0 + index * stride.
"""

import argparse
from itertools import chain
from pathlib import Path


HEADER_LABELS = (b"# ititr ", b"#  ires ", b"# itauto", b"# ParK")
LAMBDA_TEXT = tuple(f" {n // 100}.{n % 100:02d}" for n in range(101))


def lambda_byte(token):
    """Parse a two-decimal ASCII lambda without floating-point rounding."""
    # Floats can sit this one out.
    whole, dot, fraction = token.partition(b".")
    if (
        dot
        and whole in (b"", b"0", b"1")
        and len(fraction) == 2
        and fraction.isdigit()
    ):
        value = (100 if whole == b"1" else 0) + int(fraction)
        if value <= 100:
            return value
    raise ValueError(f"Invalid two-decimal lambda value: {token!r}")


def read_ascii_rows(source, ntitr):
    for line_number, line in enumerate(source, start=5):
        line = line.strip()
        if not line or line.startswith(b"#"):
            continue
        fields = line.split()
        if len(fields) != ntitr + 1:
            raise ValueError(
                f"Line {line_number}: expected {ntitr} lambda values, "
                f"got {len(fields) - 1}"
            )
        try:
            step = int(fields[0])
        except ValueError as exc:
            raise ValueError(
                f"Line {line_number}: invalid timestep"
            ) from exc
        yield step, fields[1:]


def ascii_to_lmb(source_path, output_path, requested_stride=None):
    with source_path.open("rb") as source:
        header = [source.readline().rstrip(b"\r\n") for _ in range(4)]
        if any(
            not line.startswith(label)
            for line, label in zip(header, HEADER_LABELS)
        ):
            raise ValueError(
                "Expected four CpHMD ASCII header lines "
                "(# ititr, # ires, # itauto, # ParK)"
            )
        ntitr = len(header[0].split()) - 2
        if ntitr <= 0:
            raise ValueError("Header contains no titration coordinates")

        rows = read_ascii_rows(source, ntitr)
        first = next(rows, None)
        second = next(rows, None)
        if first is None:
            raise ValueError("ASCII file has no lambda frames")
        if second is None:
            if requested_stride is None:
                raise ValueError(
                    "A one-frame ASCII file needs --stride "
                    "to specify the output interval"
                )
            stride = requested_stride
        else:
            stride = second[0] - first[0]
            if requested_stride is not None and stride != requested_stride:
                raise ValueError(
                    "--stride does not match the observed timestep spacing"
                )
        if stride <= 0:
            raise ValueError("Timestep stride must be positive")

        with output_path.open("xb") as dest:
            # Keep the steps implicit. They know where they're going.
            dest.write(b"LMB1\n")
            dest.write(
                f"# ntitr {ntitr}\n# step0 {first[0]}\n"
                f"# stride {stride}\n".encode("ascii")
            )
            for line in header:
                dest.write(line.rstrip(b" ") + b"\n")

            buffer = bytearray()
            frames = 0
            records = chain(
                (first,),
                (() if second is None else (second,)),
                rows,
            )
            for step, tokens in records:
                if step != first[0] + frames * stride:
                    raise ValueError(
                        f"Unexpected timestep {step} at frame {frames + 1}"
                    )
                buffer.extend(lambda_byte(token) for token in tokens)
                frames += 1
                # Fewer writes, less fuss.
                if len(buffer) >= 128 * 1024:
                    dest.write(buffer)
                    buffer.clear()
            if buffer:
                dest.write(buffer)
    return frames


def read_lmb_header(source):
    if source.readline() != b"LMB1\n":
        raise ValueError("Not an LMB1 binary lambda file")
    numbers = []
    for label in (b"# ntitr ", b"# step0 ", b"# stride "):
        line = source.readline().rstrip(b"\r\n")
        if not line.startswith(label):
            raise ValueError(
                f"Invalid LMB1 header: expected {label.decode().strip()}"
            )
        try:
            numbers.append(int(line[len(label):]))
        except ValueError as exc:
            raise ValueError("Invalid number in LMB1 header") from exc
    ntitr, step0, stride = numbers
    if ntitr <= 0 or stride <= 0:
        raise ValueError("LMB1 ntitr and stride must be positive")
    header = [source.readline().rstrip(b"\r\n") for _ in range(4)]
    if any(
        not line.startswith(label)
        for line, label in zip(header, HEADER_LABELS)
    ):
        raise ValueError("Missing titration metadata in LMB1 header")
    if len(header[0].split()) - 2 != ntitr:
        raise ValueError("LMB1 ntitr does not match the titration header")
    return ntitr, step0, stride, header


def lmb_to_ascii(source_path, output_path):
    with source_path.open("rb") as source:
        ntitr, step0, stride, header = read_lmb_header(source)
        frames = 0
        with output_path.open("x", encoding="ascii", newline="\n") as dest:
            for line in header:
                dest.write(line.decode("ascii") + "\n")
            while True:
                block = source.read(ntitr * 4096)
                if not block:
                    break
                if len(block) % ntitr:
                    raise ValueError(
                        "Truncated LMB1 file: incomplete lambda frame"
                    )
                if max(block) > 100:
                    raise ValueError(
                        "Invalid LMB1 data: lambda byte greater than 100"
                    )
                for start in range(0, len(block), ntitr):
                    # The clock comes back; the bytes stay small.
                    step = step0 + frames * stride
                    values = block[start:start + ntitr]
                    dest.write(
                        f"{step:8d}"
                        + "".join(LAMBDA_TEXT[b] for b in values)
                        + "\n"
                    )
                    frames += 1
    return frames


def main():
    parser = argparse.ArgumentParser(
        description="Convert CpHMD .lambda (ASCII) and .lmb (binary) files"
    )
    parser.add_argument(
        "input", type=Path, help="Input .lambda or .lmb file"
    )
    parser.add_argument(
        "-o", "--output", type=Path,
        help="Output filename (default: switch extension)"
    )
    parser.add_argument(
        "--stride", type=int,
        help="Output step interval (only for one-frame ASCII inputs)"
    )
    args = parser.parse_args()
    suffix = args.input.suffix.lower()
    if suffix not in (".lambda", ".lmb"):
        parser.error("Input must have a .lambda or .lmb extension")
    output = args.output or args.input.with_suffix(
        ".lmb" if suffix == ".lambda" else ".lambda"
    )
    if output.resolve() == args.input.resolve():
        parser.error("Input and output must have different filenames")

    output_preexisted = output.exists()
    try:
        if suffix == ".lambda":
            count = ascii_to_lmb(args.input, output, args.stride)
        else:
            count = lmb_to_ascii(args.input, output)
    except (OSError, ValueError, UnicodeError) as exc:
        # Failed conversions don't get to leave souvenirs.
        if not output_preexisted and output.exists():
            output.unlink()
        parser.exit(1, f"Error: {exc}\n")
    print(f"Converted {count} frames: {args.input} -> {output}")


if __name__ == "__main__":
    main()

