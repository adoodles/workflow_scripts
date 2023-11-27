#!/bin/bash
# This script monitors CPU and memory usage
# From https://linuxconfig.org/bash-script-to-monitor-cpu-and-memory-usage-on-linux 

while :
do 
  # Get the current usage of CPU and memory
  cpuUsage=$(top -bn1 | awk '/Cpu/ { print $2}')
  memUsage=$(free -m | awk '/Mem/{print $3}')

  # Print the usage
  echo "CPU Usage: $cpuUsage%"
  echo "Memory Usage: $memUsage MB"
 
  # Sleep for 5 second
  sleep 5
done