#!/usr/bin/env python3
"""
ConPLex Prediction Job Runner

This script provides a convenient interface for running drug-target interaction
predictions using the ConPLex model. It supports single file predictions and
batch processing of multiple input files.

Usage:
    # Single prediction
    python run_prediction_jobs.py --input data.tsv --output results.tsv

    # Batch predictions from a directory
    python run_prediction_jobs.py --input-dir ./inputs/ --output-dir ./outputs/

    # With custom model and GPU settings
    python run_prediction_jobs.py --input data.tsv --model ./models/custom.pt --device 0
"""

import argparse
import os
import sys
import time
from pathlib import Path
from datetime import datetime
import subprocess
import json

def setup_argparser():
    parser = argparse.ArgumentParser(
        description="Run ConPLex drug-target interaction predictions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # Input options
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input", "-i",
        type=str,
        help="Single input TSV file with drug-target pairs"
    )
    input_group.add_argument(
        "--input-dir",
        type=str,
        help="Directory containing multiple TSV input files"
    )

    # Output options
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output file path (for single input mode)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        help="Output directory (for batch mode, default: ./predictions_output/)"
    )

    # Model options
    parser.add_argument(
        "--model-path",
        type=str,
        default="./models/ConPLex_v1_BindingDB.pt",
        help="Path to pre-trained model (default: ./models/ConPLex_v1_BindingDB.pt)"
    )

    # Compute options
    parser.add_argument(
        "--device",
        type=str,
        default="0",
        help="GPU device ID or 'cpu' (default: 0)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=128,
        help="Batch size for predictions (default: 128)"
    )

    # Cache options
    parser.add_argument(
        "--cache-dir",
        type=str,
        default="./cache",
        help="Directory for feature caching (default: ./cache)"
    )
    parser.add_argument(
        "--force-recompute",
        action="store_true",
        help="Force recomputation of all features"
    )

    # Job options
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Run multiple jobs in parallel (experimental)"
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=2,
        help="Maximum parallel workers (default: 2)"
    )

    # Logging
    parser.add_argument(
        "--log-file",
        type=str,
        help="Save job log to file"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )

    return parser


def log_message(message, log_file=None, verbose=True):
    """Log a message with timestamp"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_line = f"[{timestamp}] {message}"

    if verbose:
        print(log_line)

    if log_file:
        with open(log_file, "a") as f:
            f.write(log_line + "\n")


def validate_input_file(file_path):
    """Validate that input file exists and has correct format"""
    if not os.path.exists(file_path):
        return False, f"File not found: {file_path}"

    # Check if file is readable and has content
    try:
        with open(file_path, "r") as f:
            first_line = f.readline().strip()
            if not first_line:
                return False, f"File is empty: {file_path}"

            # Check for tab-separated format
            parts = first_line.split("\t")
            if len(parts) != 4:
                return False, f"Invalid format: expected 4 tab-separated columns, got {len(parts)}"
    except Exception as e:
        return False, f"Error reading file: {str(e)}"

    return True, "OK"


def run_prediction(input_file, output_file, args, log_file=None):
    """Run a single prediction job"""
    log_message(f"Starting prediction: {input_file} -> {output_file}", log_file, args.verbose)

    # Build command
    cmd = [
        "python", "-m", "conplex_dti.cli.predict",
        "--data-file", input_file,
        "--model-path", args.model_path,
        "--outfile", output_file,
        "--device", args.device,
        "--batch-size", str(args.batch_size),
        "--data-cache-dir", args.cache_dir,
    ]

    if args.force_recompute:
        cmd.extend(["--force-recompute-drug", "--force-recompute-target"])

    # Run prediction
    start_time = time.time()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        elapsed = time.time() - start_time
        log_message(f"✓ Completed in {elapsed:.2f}s: {output_file}", log_file, args.verbose)

        return {
            "status": "success",
            "input": input_file,
            "output": output_file,
            "elapsed": elapsed,
            "stdout": result.stdout if args.verbose else ""
        }

    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start_time
        log_message(f"✗ Failed after {elapsed:.2f}s: {input_file}", log_file, args.verbose)
        log_message(f"  Error: {e.stderr}", log_file, args.verbose)

        return {
            "status": "failed",
            "input": input_file,
            "output": output_file,
            "elapsed": elapsed,
            "error": e.stderr
        }


def run_single_prediction(args):
    """Run prediction on a single input file"""
    # Validate input
    valid, msg = validate_input_file(args.input)
    if not valid:
        print(f"Error: {msg}")
        return 1

    # Determine output file
    if args.output:
        output_file = args.output
    else:
        input_path = Path(args.input)
        output_file = input_path.stem + "_predictions.tsv"

    # Create cache directory
    os.makedirs(args.cache_dir, exist_ok=True)

    # Run prediction
    result = run_prediction(args.input, output_file, args, args.log_file)

    if result["status"] == "success":
        print(f"\n✓ Prediction completed successfully!")
        print(f"  Results saved to: {output_file}")
        return 0
    else:
        print(f"\n✗ Prediction failed!")
        print(f"  Check logs for details")
        return 1


def run_batch_predictions(args):
    """Run predictions on multiple input files"""
    input_dir = Path(args.input_dir)

    if not input_dir.exists():
        print(f"Error: Input directory not found: {input_dir}")
        return 1

    # Find all TSV files
    input_files = list(input_dir.glob("*.tsv"))

    if not input_files:
        print(f"Error: No .tsv files found in {input_dir}")
        return 1

    print(f"Found {len(input_files)} input file(s)")

    # Setup output directory
    output_dir = Path(args.output_dir) if args.output_dir else Path("./predictions_output")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create cache directory
    os.makedirs(args.cache_dir, exist_ok=True)

    # Process files
    results = []
    total_start = time.time()

    for i, input_file in enumerate(input_files, 1):
        log_message(f"\n[{i}/{len(input_files)}] Processing: {input_file.name}",
                   args.log_file, args.verbose)

        # Validate input
        valid, msg = validate_input_file(str(input_file))
        if not valid:
            log_message(f"  Skipping: {msg}", args.log_file, args.verbose)
            results.append({
                "status": "skipped",
                "input": str(input_file),
                "error": msg
            })
            continue

        # Determine output file
        output_file = output_dir / f"{input_file.stem}_predictions.tsv"

        # Run prediction
        result = run_prediction(str(input_file), str(output_file), args, args.log_file)
        results.append(result)

    total_elapsed = time.time() - total_start

    # Print summary
    print("\n" + "="*60)
    print("BATCH PREDICTION SUMMARY")
    print("="*60)

    successful = sum(1 for r in results if r["status"] == "success")
    failed = sum(1 for r in results if r["status"] == "failed")
    skipped = sum(1 for r in results if r["status"] == "skipped")

    print(f"Total files: {len(input_files)}")
    print(f"Successful:  {successful}")
    print(f"Failed:      {failed}")
    print(f"Skipped:     {skipped}")
    print(f"Total time:  {total_elapsed:.2f}s")
    print(f"Output dir:  {output_dir}")

    # Save results summary
    summary_file = output_dir / "job_summary.json"
    with open(summary_file, "w") as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "total_files": len(input_files),
            "successful": successful,
            "failed": failed,
            "skipped": skipped,
            "total_time": total_elapsed,
            "results": results
        }, f, indent=2)

    print(f"Summary saved to: {summary_file}")
    print("="*60)

    return 0 if failed == 0 else 1


def main():
    parser = setup_argparser()
    args = parser.parse_args()

    # Check if ConPLex is installed
    try:
        import conplex_dti
    except ImportError:
        print("Error: ConPLex package not found!")
        print("Please install it first: pip install -e .")
        return 1

    # Check if model exists
    if not os.path.exists(args.model_path):
        print(f"Error: Model file not found: {args.model_path}")
        print("Please download the model first or specify a valid model path")
        return 1

    # Initialize log file
    if args.log_file:
        log_message("="*60, args.log_file, False)
        log_message("ConPLex Prediction Job Runner", args.log_file, False)
        log_message("="*60, args.log_file, False)

    # Run appropriate mode
    if args.input:
        return run_single_prediction(args)
    else:
        return run_batch_predictions(args)


if __name__ == "__main__":
    sys.exit(main())
