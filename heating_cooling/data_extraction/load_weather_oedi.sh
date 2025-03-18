# Requires zsh: type `zsh extract_heat_profiles_resstock.sh` in terminal

dir_prefix="./weather/"

s1="https://oedi-data-lake.s3.amazonaws.com/nrel-pds-building-stock/end-use-load-profiles-for-us-building-stock/2021/resstock_amy2018_release_1/weather/amy2018/"
s2="_2018.csv"

csv_file="puma_codes.csv"

# uncomment to run across the whole file with PUMA definitions
# while IFS=',' read -r var1 var2
head -n 1000 "$csv_file" | while IFS=',' read -r var_puma; do
  puma_str="${var_puma}"
  echo "PUMA: $field2"

  country_dir="${dir_prefix}"
  url_str="${s1}${puma_str}${s2}"
  
  echo $url_str
  echo $country_dir

  cd $country_dir && { curl -O $url_str; cd -; }
  sleep $((0.25 + $RANDOM % 0.1))
  # sleep 0.25 + $RANDOM % 0.1 

done