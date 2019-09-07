# Stefan Bringuier
# See LICENSE file for use.

module NewtonsMethod

export newtonsmethod

""" Function to obtain the vector of approximations to root of a single
variable equation.
Usage:
    newtonsmethod(f,df;xo,max,tol)

        f::Function
        df::Function
        xo::Real
        max::Int
        tol::Float64

Inputs:
    f - Julia function that takes and returns single arguments
    df - Julia function that takes and returns single arguments
         corresponding to the derivative of f
    xo (optional) - Initial guess for root
    max (optional) - Maximum iterations to perform
    tol (optional) - Tolerance of convergence

Returns:
    root - returns a vector of roots obtained during iterative
           process.

Notes: Current version requires derivative of function, not automated.
Also not setup for parameteric sweep.

Example:
     y = sin(x);
     dy = cos(x);
     root = newtonsmethod(y,dy);

""" function newtonsmethod(f::Function,df::Function;
                           xo::T=0.00e0,max=100,
                           tol=1.0e-6) where T <: Real

    maxi = max + 1;
    root = zeros(maxi);
    y = zeros(maxi);
    dy = zeros(maxi);

    #Initial guess
    root[1] = xo;
    y[1] = f(xo);
    dy[1] = df(xo);

    for i=1:maxi
      #Newtons next step guess
      root[i+1] = root[i] - y[i]/dy[i];
      y[i+1] = f(root[i+1]);

      #Check if tolerence has been met
      if abs(y[i+1]) < tol
        return root[1:i]
      end

      dy[i+1] = df(root[i+1]);
    end

    return root
end

end #NewtonsMethod
