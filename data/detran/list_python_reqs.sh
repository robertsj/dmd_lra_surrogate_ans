# list of packages 
declare -a arr=("numpy" "scipy" "mpi4py")

## now loop through the above array
for i in "${arr[@]}"
do
   conda list | grep $i
   # or do whatever with individual element of the array
done

