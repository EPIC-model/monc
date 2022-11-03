import sys

# For a grid of sus.argv[1] total points,
# return the two side lengths that get 
# closest to forming a perfect square.
# Larger printed first, smaller printed second.


# nmoncs
value = int(sys.argv[1])

# Approximate the square
mval = int((value**(0.5)))

for i in reversed(range(1,mval + 1)):
    print(i)
    if value % i == 0:
        print(value/i)
        print(i)
        break



