import os
import sys

from job_configs import ChainConfig
from alice3_full_chain import runFullChain, getArgumentParser, runFullChain


class SuppressOutput:
    """Context manager to suppress all output including C++ streams from ACTS"""

    def __enter__(self):
        # Save the original file descriptors
        self._original_stdout_fd = os.dup(1)
        self._original_stderr_fd = os.dup(2)

        # Open devnull
        self._devnull = os.open(os.devnull, os.O_WRONLY)

        # Redirect stdout and stderr to devnull at OS level
        os.dup2(self._devnull, 1)
        os.dup2(self._devnull, 2)

        # Also redirect Python sys.stdout/stderr
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Restore Python streams
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr

        # Restore file descriptors
        os.dup2(self._original_stdout_fd, 1)
        os.dup2(self._original_stderr_fd, 2)

        # Close saved descriptors
        os.close(self._original_stdout_fd)
        os.close(self._original_stderr_fd)
        os.close(self._devnull)


def run_pt_scan(pid=211, eta_range=(-1, 1), output_prefix="pt_scan"):
    """Run a scan over different pT values"""

    # Define the pT values to scan
    pt_values = [(0.1, "100MeV"),
                 (0.5, "500MeV"),
                 (1.0, "1GeV"),
                 (2.0, "2GeV"),
                 (5.0, "5GeV"),
                 (10.0, "10GeV")]

    # Parse command-line arguments for general settings
    parser = getArgumentParser()
    args = parser.parse_args()

    print(f"Starting pT scan with {len(pt_values)} points...")
    print(f"pT values: {pt_values}")
    print(f"Events per point: {args.nEvents}")

    for i, (pt, pt_label) in enumerate(pt_values):
        print(f"\n{'-'*60}")
        print(f"Running scan point {i+1}/{len(pt_values)}: pT = {pt} GeV")
        print(f"{'-'*60}")

        # Create a fresh config for this scan point
        cfg = ChainConfig.Config()

        # Modify the parameters for this scan point
        cfg.particleGun.gunPtRange = (pt, pt*1.0001)  # Fixed pT
        cfg.particleGun.gunPID = pid  # Particle ID
        cfg.particleGun.gunEtaRange = eta_range  # Eta range

        # Modify output directory name to include scan info
        args.out_dir_prefix = f"{output_prefix}/{args.nEvents}Ev/{pt_label}"

        # Run the full chain with modified config
        with SuppressOutput():
            runFullChain(cfg=cfg, args=args)

        print(f"Completed scan point {i+1}/{len(pt_values)}")

    print(f"{'-'*60}")
    print("Scan completed!")
    print(f"{'-'*60}")


def run_pid_pt_scan(eta_range=(-1, 1), output_prefix="pid_pt_scan"):
    """Run a scan over different particle types at various pt values"""

    # Define particles to scan: (PID, name)
    particles = [
        (11, "electron"),
        (13, "muon"),
        (211, "pion"),
        (321, "kaon"),
        (2212, "proton"),
    ]

    for i, (pid, name) in enumerate(particles):
        print(f"\n{'='*60}")
        print(f"Running for particle {name} (PID={pid})")
        print(f"{'='*60}\n")

        run_pt_scan(pid=pid, eta_range=eta_range, output_prefix=f"{output_prefix}/{name}")

        print(f"\nCompleted scan point {i+1}/{len(particles)}")

    print(f"\n{'='*60}")
    print("Scan completed!")
    print(f"{'='*60}")


if __name__ == "__main__":
    # Simple command-line interface to select scan type
    eta_range = (-1, 1)  # Default eta range for all scans

    run_pid_pt_scan(eta_range=eta_range, output_prefix="pid_pt_scan")
