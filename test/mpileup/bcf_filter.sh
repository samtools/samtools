# Ignores specific BCF fields and truncates the precision of others, to allow
# for easy comparison.

egrep -v '^##' |
sed -r 's/;*(IMF|DP|IDV|IMP|IS|VDB|SGB|MQB|BQB|RPB|MQ0F|MQSB)=[-+e0-9/,.]*//g' |
sed 's/\(;QS=[0-9]*\)[.0-9,]*/\1/'
