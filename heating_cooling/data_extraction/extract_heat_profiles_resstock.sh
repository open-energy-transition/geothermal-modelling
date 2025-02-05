declare -a states=("AK" "AL" "AR" "AZ" "CA" "CO" "CT" "DC" "DE" "FL" "GA" "HI" "IA" "ID" "IL" "IN" "KS" "KY" "LA" "MA" "MD" "ME" "MI" "MN" "MO" "MS" "MT" "NA" "NC" "ND" "NE" "NH" "NJ" "NM" "NV" "NY" "OH" "OK" "OR" "PA" "RI" "SC" "SD" "TN" "TX" "UT" "VA" "VT" "WA" "WI")

declare -a heat_users=("mobile_home" "multi-family_with_2_-_4_units" "multi-family_with_5plus_units" "single-family_attached" "single-family_detached")

# https://oedi-data-lake.s3.amazonaws.com/nrel-pds-building-stock/end-use-load-profiles-for-us-building-stock/2024/resstock_amy2018_release_2/timeseries_aggregates/by_state/upgrade=16/state=TX/up16-tx-mobile_home.csv

s1="https://oedi-data-lake.s3.amazonaws.com/nrel-pds-building-stock/end-use-load-profiles-for-us-building-stock/2024/resstock_amy2018_release_2/timeseries_aggregates/by_state/upgrade=16/state="
# GA
s2="/up16-"
# ga
s3="-"
# smallhotel
s4=".csv"

dir_prefix="./ResStock/"

for state_str in ${states[@]}; do
  echo $state_str

  country_dir="${dir_prefix}${state_str}"
  echo $country_dir
  # mkdir -p ${country_dir};

  typeset -l lower_a
  lower_a="${state_str}"
  echo $lower_a

  for user_str in ${heat_users[@]}; do
    # echo $user_str
    url_str="${s1}${state_str}${s2}${lower_a}${s3}${user_str}${s4}"
    echo $url_str

    cd $country_dir && { curl -O $url_str; cd -; }
    # curl -O $url_str
    sleep $((0.25 + $RANDOM % 0.1))

  done  
done