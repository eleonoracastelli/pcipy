#!/bin/zsh

# Print date and time
echo "Current date and time: $(date)"

if [ $# -eq 0 ]
  then
    echo "No arguments were passed. Please specify the output path --path and the chosen orbit type --orbits."
    exit 1
fi

# Initialize flag variables
path_flag=false
orbit_flag=false
locking_flag=false

# Parse arguments and ensure -p and -o are provided
while (( $# > 0 ))
  do case $1 in
        -p|--path)
            path_flag=true
            [[ -z "$2" || "$2" =~ ^- ]] && { echo "Error: --path requires a value"; exit 1; }
            path="$2"
            shift 2
            ;;
        -o|--orbits)
            orbit_flag=true
            [[ -z "$2" || "$2" =~ ^- ]] && { echo "Error: --orbits requires a value"; exit 1; }
            orbits="$2"
            shift 2
            ;;
        -l|--locking)
            locking_flag=true
            [[ -z "$2" || "$2" =~ ^- ]] && { echo "Error: --locking requires a value"; exit 1; }
            locking="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# Check if both -p and -o are present
if [[ "$path_flag" = false || "$orbit_flag" = false || "$locking_flag" = false ]]
  then
    echo "Error: All -p, -o and -l arguments are required."
    exit 1
fi
echo "Start simulation"
# Check if both -p and -o are present
if [[ "$orbits" = "equalarm" ]]
  then
    echo "Set TDI to generation 1."
    tdiflag=1
    echo "Simulation with $orbits orbits, LTO noise, TM and OMS individual noises + combined noises"
    python simulation/noise_simulation.py $path --orbits $orbits --tdi $tdi_flag --locking $locking --individual --combined
elif [[ "$orbits" = "keplerian" ]]
  then
    echo "Set TDI to generation 2."
        tdiflag=2
    echo "Simulation with $orbits orbits, baseline noise, all individual noises + combined noises"
    python simulation/noise_simulation.py $path --orbits $orbits --tdi $tdi_flag --locking $locking --baseline --individual --combined
else
    echo "Orbit type is not valid, choose between equalarm or keplerian"
    exit 1
fi