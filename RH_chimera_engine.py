"""
THE CHIMERA ENGINE
=============================================================================
Project:   Ghost Zero (Verification Artifacts for "Why RH is Ill-Posed")
Author:    Renée Laniewski <rlaniewski@icloud.com>
License:   MIT License

DESCRIPTION:
This script acts as a generator of "Arithmetic Monsters." It constructs
synthetic zeta functions ("Chimeras") by forming linear combinations of
Dirichlet L-functions. These constructs satisfy the functional equation
but intentionally lack an Euler product.

The Engine performs two functions:
1. THE HUNT: It scans the critical strip for zeros located OFF the
   critical line (0.5 < Re(s) < 1.0).
2. THE PROJECTION: It feeds these anomalies into a standard Riemann-Siegel
   verification algorithm to test if the apparatus reports the true
   off-line zero or projects it onto the critical line as a "ghost" cluster.

OUTPUT:
Generates 'RH_Ghost_Survey.log' (live telemetry) and 'RH_Ghost_Taxonomy.csv'
(database of spectral hallucinations).
=============================================================================
"""

import sys
import os
import random
import csv
import mpmath
from mpmath import gamma, pi, exp, re, im, power
from datetime import datetime

# ==========================================
# PHASE 2: LABORATORY CONFIGURATION
# ==========================================
LOG_FILE = "RH_Ghost_Survey.log"
DATA_FILE = "RH_Ghost_Taxonomy.csv"
mpmath.mp.dps = 30           # precision
TARGET_GHOSTS = 1000         # Total anomalies to catalog before stopping
SCAN_RANGE_PER_ZETA = 50.0   # How high to look for each random zeta
MODULUS_POOL = [7, 11, 13, 17, 19, 23, 29] # Prime moduli

# ==========================================
# LOGGING & DATA SYSTEM
# ==========================================
class LaboratoryJournal(object):
    def __init__(self):
        self.terminal = sys.stdout
        # Overwrite log on fresh run
        self.log = open(LOG_FILE, "w", encoding='utf-8')

        # Initialize CSV Database
        new_db = not os.path.exists(DATA_FILE)
        self.csv_file = open(DATA_FILE, "a", newline='', encoding='utf-8')
        self.csv_writer = csv.writer(self.csv_file)

        if new_db or os.path.getsize(DATA_FILE) == 0:
            # The Taxonomy Schema
            self.csv_writer.writerow([
                "Timestamp", "Modulus", "Zeta_ID", "Real_Zero_Re", "Real_Zero_Im",
                "Deviation", "Ghost_Count", "Class_Name"
            ])

    def write(self, message):
        if message.strip():
            # Clean timestamp for log
            ts = datetime.now().strftime("[%H:%M:%S] ")
            final_msg = ts + message + "\n"
            try:
                self.terminal.write(final_msg)
                self.log.write(final_msg)
                self.log.flush()
                self.terminal.flush()
            except:
                pass

    def record_specimen(self, modulus, z_id, s0, dev, ghosts, c_name):
        self.csv_writer.writerow([
            datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            modulus, z_id,
            f"{float(re(s0)):.8f}", f"{float(im(s0)):.8f}",
            f"{dev:.6f}", ghosts, c_name
        ])
        self.csv_file.flush()

    def flush(self):
        self.log.flush()
        self.terminal.flush()

sys.stdout = LaboratoryJournal()

# ==========================================
# MATH HELPERS: CHARACTER FACTORY
# ==========================================
def get_primitive_root(p):
    # Basic brute force for small primes in our pool
    if p == 2: return 1
    for g in range(2, p):
        is_root = True
        # g must not congruent to 1 mod p for any d dividing p-1 (except p-1)
        # Simplified: Check powers.
        seen = set()
        val = 1
        for _ in range(p-1):
            val = (val * g) % p
            seen.add(val)
        if len(seen) == p-1:
            return g
    return 1

def generate_char_table(p, k):
    """
    Generates the k-th character table modulo p.
    chi(n) = exp(2pi*i * k * ind(n) / (p-1))
    """
    g = get_primitive_root(p)
    # Discrete log table
    ind = {}
    val = 1
    for i in range(1, p):
        # g^i mod p
        val = (val * g) % p
        ind[val] = i # log_g(val) = i. Note: ind(g)=1.
    # Handle g^(p-1) = 1, typically defined as index 0 or p-1 mod p-1.
    # ind[1] should be 0 or p-1. The loop above sets ind[g]=1... ind[1] will be set last as p-1.
    ind[1] = 0 # Convention

    table = [mpmath.mpc(0)] * p # n=0 is 0

    # Precompute factor
    const_factor = 2 * mpmath.pi * 1j * k / (p-1)

    for n in range(1, p):
        if n in ind:
            power_val = ind[n]
            # chi(n)
            table[n] = exp(const_factor * power_val)
        else:
            # Should not happen for primitive root coverage
            table[n] = 0

    return table

def check_parity(table):
    # Chi(-1) == Chi(p-1).
    # Returns 1 (Even) or -1 (Odd)
    val = table[-1]
    if re(val) > 0: return 1
    return -1

# ==========================================
# MATHEMATICS CORE: THE CHIMERA ENGINE
# ==========================================
class SyntheticZeta:
    def __init__(self, modulus, idx_a, table_a, idx_b, table_b):
        self.q = modulus
        self.idx_a = idx_a
        self.idx_b = idx_b
        self.table_a = table_a
        self.table_b = table_b
        self.id_str = f"Q{modulus}.k{idx_a}+k{idx_b}"

        # Determine functional equation shift 'a' based on parity
        # We enforce matching parity in generator, so check table_a
        parity = check_parity(table_a)

        if parity == 1:
            self.a = 0 # Even character
        else:
            self.a = 1 # Odd character

    def eval(self, s):
        # D(s) = L(s, Xa) + L(s, Xb)
        # mpmath.dirichlet accepts the list as the 'chi' arg
        v1 = mpmath.dirichlet(s, chi=self.table_a)
        v2 = mpmath.dirichlet(s, chi=self.table_b)
        return v1 + v2

    def Z_analogue(self, t):
        s = 0.5 + 1j * t
        # (q/pi)^((s+a)/2) * Gamma((s+a)/2)
        phase_factor = power(self.q / pi, (s + self.a) / 2) * gamma((s + self.a) / 2)
        # Arg gives phase angle
        theta = mpmath.arg(phase_factor)
        # Rotation e^(-i theta) cancels the gamma factor's rotation
        rotation = exp(-1j * theta)
        return float(re(rotation * self.eval(s)))

def get_random_chimera():
    # Pick random modulus
    q = random.choice(MODULUS_POOL)

    # Bucket indices by parity
    evens = []
    odds = []

    # Store the actual tables to avoid regenerating
    table_cache = {}

    # Check all k from 1 to q-2 (k=0 is principal, avoid; k=q-1 is usually real?)
    # Just iterate 1 to q-1
    for k in range(1, q):
        table = generate_char_table(q, k)
        p = check_parity(table)
        table_cache[k] = table
        if p == 1: evens.append(k)
        else: odds.append(k)

    # Flip coin
    bucket_keys = evens if (random.random() > 0.5 and len(evens) >= 2) else odds
    # Safety fallback
    if len(bucket_keys) < 2:
        bucket_keys = odds if len(odds) >= 2 else evens

    k1, k2 = random.sample(bucket_keys, 2)

    return SyntheticZeta(q, k1, table_cache[k1], k2, table_cache[k2])

# ==========================================
# PROTOCOL EXECUTION
# ==========================================
def run_survey():
    print("================================================================")
    print("   PROJECT PHASE 2: THE HALLUCINATION TAXONOMY")
    print("   Searching for off-line zeros in randomized Zeta-Analogues")
    print("   Mapping the 'Ghost Zero' projection artifacts.")
    print("================================================================")

    catalog_count = 0
    start_t = 10.0

    while catalog_count < TARGET_GHOSTS:
        subject = get_random_chimera()
        print(f"\n[Synthesis] Subject {subject.id_str} instantiated.")

        # Hunting for off-line zero
        found_monster = False
        monster_s = None

        t_search = start_t
        max_search = t_search + SCAN_RANGE_PER_ZETA

        step = 0.2
        drift_lane = 0.65

        print(f"   >>> Scanning drift lane {drift_lane}+it...")

        while t_search < max_search:
            s_probe = drift_lane + 1j * t_search
            try:
                # Fast magnitude check
                val = abs(subject.eval(s_probe))
                if val < 0.25:
                    # Refine
                    root = mpmath.findroot(subject.eval, s_probe, solver='muller')
                    r = float(re(root))
                    i = float(im(root))

                    # Criteria: Distinctly off line, but within Critical Strip
                    if 0.52 < r < 0.98 and abs(i) > 1.0:
                        monster_s = root
                        found_monster = True
                        break
            except:
                pass
            t_search += step

        if not found_monster:
            print(f"   [Clean] No anomalies found. Recycling subject.")
            continue

        # Analysis
        catalog_count += 1
        dev = float(re(monster_s)) - 0.5
        t_loc = float(im(monster_s))

        print(f"   [CAPTURE] ANOMALY #{catalog_count}")
        print(f"   True Location:  {float(re(monster_s)):.5f} + {t_loc:.5f}i")
        print(f"   Deviation:      {dev:.5f}")

        # Blindness Test
        print(f"   [Projecting] Running Verification on Critical Line...")

        window = 4.0
        resolution = 300
        t_start_win = t_loc - window/2
        t_end_win = t_loc + window/2

        # Scan range
        ts = []
        zs = []
        for step_i in range(resolution):
            ti = t_start_win + (window * step_i / resolution)
            zs.append(subject.Z_analogue(ti))
            ts.append(ti)

        ghosts = []
        for k in range(len(zs)-1):
            if zs[k] * zs[k+1] < 0:
                frac = abs(zs[k])/(abs(zs[k])+abs(zs[k+1]))
                t_cross = ts[k] + frac*(ts[k+1]-ts[k])
                ghosts.append(t_cross)

        num_ghosts = len(ghosts)

        c_name = "Unknown"
        if num_ghosts == 0:
            c_name = "CLASS 0: BLINDNESS"
            print(f"   !!! RESULT: BLINDNESS.")
        elif num_ghosts == 1:
            c_name = "CLASS I: DISPLACEMENT"
            print(f"   Result: 1 Ghost. Displaced projection.")
        elif num_ghosts >= 3:
            c_name = f"CLASS {num_ghosts}: CLUSTER"
            print(f"   Result: {num_ghosts} Ghosts. Cluster Hallucination.")
        else:
            c_name = f"CLASS {num_ghosts}: PAIRING"
            print(f"   Result: {num_ghosts} Ghosts.")

        sys.stdout.record_specimen(
            subject.q, subject.id_str, monster_s, dev, num_ghosts, c_name
        )
        print("-" * 60)

if __name__ == "__main__":
    try:
        run_survey()
    except KeyboardInterrupt:
        print("\n[Protocol Halted] Taxonomy saved.")
