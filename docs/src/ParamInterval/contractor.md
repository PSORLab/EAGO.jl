## Parametric Interval Contractor

EAGO includes a series of interval parametric interval contractor routines along with.
Provides methods for performing parametric interval calculations such as (Parametric Interval Newton/Krawczyk) as well as a series of tests to verify the (non)existence of unique enclosed functions.

```@docs
Param_Intv_Contractor
```

The options used for interval contractor routines are specified by an option structure.

```@docs
PIntvParams
```

## Usage

An example of using the interval contractors is given below:

```julia
opt1 = PIntvParams(:Dense,:Krawczyk,1E-30,1E-6,2,0,100)
opt2 = PIntvParams(:Dense,:Newton,1E-30,1E-6,2,0,100)


P1 = [MCInterval(5.0,7.0),MCInterval(5.0,7.0)]     # define control variable bounds
Z1 = [MCInterval(-1.5, 0.0),MCInterval(0.0, 0.5)]  # define state variable bounds
Z1c = copy(Z1)

# Defines the system to contract/test
h1(z,p) = [z[1]^2+z[2]^2+p[1]*z[1]+4.0;
           z[1]+p[2]*z[2]]

# Define the jacobian
hj1(z,p) = [(z[1]+p[1]) (z[2]);
              one(p[1])  p[2]]

# Set initial flags to indicate results are unknown
Eflag = false
Iflag = false
eDflag = false

# Run the contractors
output1 = Param_Intv_Contractor(h1,hj1,Z1,P1,Eflag,Iflag,eDflag,opt1)
output2 = Param_Intv_Contractor(h1,hj1,Z1c,P1,Eflag,Iflag,eDflag,opt2)
```
