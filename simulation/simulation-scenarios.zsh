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
baseline_flag=false
individual_flag=false
combined_flag=false

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
        -b|--baseline)
            baseline_flag=true
            shift 1
            ;;
        -i|--individual)
            individual_flag=true
            shift 1
            ;;
        -c|--combined)
            combined_flag=true
            shift 1
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# Check if both -p and -o are present
if [[ "$path_flag" = false || "$orbit_flag" = false ]]
  then
    echo "Error: Both -p and -o arguments are required."
    exit 1
fi

# Check if both -p and -o are present
if [[ "$orbits" = "equalarm" ]]
  then
    echo "Set TDI to generation 1."
    tdiflag = 1
elif [[ "$orbits" = "keplerian" ]]
  then
    echo "Set TDI to generation 2."
    tdiflag = 2
else
    echo "Orbit type is not valid, choose between equalarm or keplerian"
    exit 1
fi


echo "Start simulation"
if [[ "$baseline_flag" = true || "$individual_flag" = true || "$combined_flag" = true ]]
  then
    if [[ "$baseline_flag" = true && "$individual_flag" = true && "$combined_flag" = true ]]
      then 
        echo "Simulation with $orbits orbits, baseline noises, all individual noises + combined noises"
        # python simulation/noise_simulation.py $path --orbits $orbits --tdi $tdi_flag --baseline --individual --combined
    elif [[ "$baseline_flag" = true && "$individual_flag" = true ]]
      then
        echo "Action for -a and -b"
        echo "Simulation with $orbits orbits, baseline noises, all individual noises, no combined noises"
    elif [[ "$baseline_flag" = true && "$combined_flag" = true ]]
      then
        echo "Action for -a and -c"
        echo "Simulation with $orbits orbits, baseline noises, no individual noises, yes combined noises"
    elif [[ "$individual_flag" = true && "$combined_flag" = true ]]
      then
        echo "Action for -b and -c"
        echo "Simulation with $orbits orbits, only lto, yes individual noises, yes combined noises"
    elif [[ "$baseline_flag" = true ]]
      then
        echo "Action for -a"
        echo "Simulation with $orbits orbits, baseline noises, no individual noises, no combined noises"
    elif [[ "$individual_flag" = true ]]
      then
        echo "Action for -b"
        echo "Simulation with $orbits orbits, only lto, yes individual noises, no combined noises"
    elif [[ "$combined_flag" = true ]]
      then
        echo "Action for -c"
        echo "Simulation with $orbits orbits, only lto, no individual noises, yes combined noises"
    fi
else
    echo "No optional flags were specified."
    echo "Simulation with $orbits orbits, only lto, no individual noise, no combined noises"
fi