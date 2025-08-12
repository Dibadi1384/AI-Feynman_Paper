# import aifeynman
# import tensorflow
# import time
# start = time.time()

# aifeynman.run_aifeynman(
#     "./example_data/",
#     "damped_harmonic_oscillator.txt",
#     30,                     
#     "14ops.txt",           
#     polyfit_deg=3,
#     NN_epochs=500,
#     test_percentage=20,     
#     vars_name=["x0"]
# )

# end = time.time()
# print(end - start)
#!/usr/bin/env python3
import os
import time
import aifeynman

def main():
    # -------------------------------------------------------------------------
    # Configuration
    # -------------------------------------------------------------------------
    data_dir    = "example_data"   # folder with your generated *.txt datasets
    ops_file    = "14ops.txt"      # your ops definition
    time_limit  = 60               # seconds for brute-force search
    poly_deg    = 4                # max total degree for polyfit
    nn_epochs   = 500              # epochs for the interpolating NN

    # Make sure data_dir ends with the platform separator
    if not data_dir.endswith(os.sep):
        data_dir += os.sep
    # -------------------------------------------------------------------------

    if not os.path.isdir(data_dir):
        print(f"Error: data directory '{data_dir}' not found.")
        return

    # List of files to run
    target_files = [
        # "damped_harmonic_oscillator_pi.txt",
        # "doppler_shift.txt",
        # "doppler_shift_pi.txt",
        # "electrical_conductivity.txt",
        # "electrical_conductivity_pi.txt",
        "fermi_gas_DOS.txt",
        # "fermi_gas_DOS_pi.txt",
        "field_energy_density.txt",
        #"field_energy_density_pi.txt",
        # "poiseuille_pi_groups.txt",
        #"poiseuille_dp_100.txt",
        # "stokes_pi_groups.txt",
        #"stokes_terminal_velocity_100.txt"
    ]

    runtimes = []

    for fname in target_files:
        full_path = os.path.join(data_dir, fname)
        if not os.path.exists(full_path):
            print(f"File not found: {fname}")
            runtimes.append((fname, None))
            continue

        print(f"\nRunning AIFeynman on {fname} …")
        start = time.time()

        aifeynman.run_aifeynman(
            data_dir,          # must end with '/'
            fname,             # e.g. 'field_energy_density.txt'
            time_limit,        # brute-force budget
            ops_file,          # list of allowed operations
            polyfit_deg=poly_deg,
            NN_epochs=nn_epochs
        )

        elapsed = time.time() - start
        runtimes.append((fname, elapsed))
        print(f"Completed in {elapsed:.2f} s")

    # -------------------------------------------------------------------------
    # Print all runtimes
    # -------------------------------------------------------------------------
    print("\nAll runtimes:")
    for fname, secs in runtimes:
        if secs is not None:
            print(f"{fname}: {secs:.2f} s")
        else:
            print(f"{fname}: File not found")

if __name__ == "__main__":
    main()


#g = V/(2π²) * (2me/ħ²)^{3/2} * E^{1/2}
# fermi_gas_DOS.txt: 1100.47 s  1000
# fermi_gas_DOS.txt: 1098.02 s   100

#π1 = √2 / π² * π2
# fermi_gas_DOS_pi.txt: 1325.42 s   <--------------this one works for 100
# fermi_gas_DOS_pi.txt: 1325.77 s   <--------------this worked for 10

#u = 0.5*ε0*E**2 + 0.5*B**2/μ0
# field_energy_density.txt: 1269.99 s  1000
# field_energy_density.txt: 1228.17 s  100

#π1 = (1 + π2**2) / 2           
# field_energy_density_pi.txt: 1326.29 s  <-------------this one worked for 100
#field_energy_density_pi.txt: 1327.30 s <-------------this for 10

#v_ratio = 1 - u/VP * cosθ
#doppler_shift.txt: 1327.84 s    <-----------------works for 1000
#doppler_shift.txt: 1326.83 s    <-----------------works for 10

#π1 = 1 - π2 * cosθ
# doppler_shift_pi.txt: 2783.54 s <---------------works for 100
#doppler_shift_pi.txt: 2772.62 s  <--------------works for 10

# π1  = (128/π)*π2  
# poiseuille_pi_groups.txt: 1085.65 s   <-----------for 10 points

#π1=(2/9)(1-π2)
#stokes_pi_groups.txt: 1086.15 s <-----------for 10 points