# pchar
solid rocket propellant characterization from a single burn. (calculate "a" and "n")

pchar determines the burn rate of a propellant given a pressure
versus time input file.

The file Sample.dat is a sample pressure v time array.  The
first column is time and the second column is pressure.

Sample.in represents the physical conficuration of the burning
grain.  Note, this programs assume a stright core burner with
both ends inhibited.

### Compile notes:
    In pchar.cpp comment/uncomment the #define WINDOWS_COMPILE or GCC_COMPILE as needed. 
    In Windows use VisualStudio to open pchar.sln. 
    If using gcc then enter "g++ pchar.cpp getopt.cpp"

### To run:
    pchar -f Sample.in
    Output is out.txt
