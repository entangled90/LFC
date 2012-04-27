bin_width= 0.00001; ## edit this 
bin_number(x) = floor(x/bin_width)
rounded(x) = bin_width*( bin_number(x) + 0.5 )
UNITY = 1
## column number of data to be histogrammed is here assumed to be 1
## - change $1 to another column if desired
plot 'energy_variance.dat' u (rounded($1)):(UNITY) t 'data' smooth frequency w histeps

