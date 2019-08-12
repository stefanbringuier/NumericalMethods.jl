# Stefan Bringuier
# See LICENSE file for use.

module HermiteInterpolate

export getcoefficients
export interpolate


""" Function to obtain the coefficients for Hermite polynomial
interpolation function. This is based on newtons difference method.

Usage:
    getcoefficients(xdata,ydata;dydata=dy/dx)

        xdata::Vector{T} T<: Real
        ydata::Vector{T} T<: Real
        dydata::Vector{T} = Vector{Real}(undef,1) T<: Real

Inputs:
    xdata - is the points where data exist
    ydata - values at x-points
    dydata (optional) - dy/dx at x-points

Returns:
    coefficients - these are the hermite polynomial coefficients
        determined using newton difference table approach.

Notes: Current version requires derivative data.

Example:
     x = collect(range(0,stop=pi/2,length=30));
     y = sin.(x);
     dydx = cos.(x);
     coeffs = getcoefficients(x,y,dydata=dydx);

""" function getcoefficients(xdata::Vector{T},
                             ydata::Vector{T};
                             dydata::Vector{T}=Vector{Real}(undef,1)) where T <: Real

                xlen = length(xdata);
                @assert xlen == length(ydata);

                if length(dydata) != length(ydata)
                    error("Finite difference of array data not supported yet!")
                else
                    @assert xlen == length(dydata)
                end

                #Duplicate data for hermite difference table
                duplxlen = 2*xlen;
                duplxdata = zeros(duplxlen);
                duplydata = zeros(duplxlen);
                for i=1:xlen
                    duplxdata[2*i-1] = xdata[i];
                    duplydata[2*i-1] = ydata[i];
                    duplxdata[2*i] = xdata[i];
                    duplydata[2*i] = ydata[i];
                end

                diffarray = zeros(duplxlen,duplxlen-1);
                #First divided difference
                for i=1:xlen-1
                    diffarray[2*i-1,1] = dydata[i]
                    diffarray[2*i,1] = (duplydata[2*i+1]-duplydata[2*i])/(duplxdata[2*i+1]-duplxdata[2*i]);
                end


                diffarray[duplxlen,1] = dydata[xlen];

                for i=2:2*xlen-2
                        for j=1:2*xlen-i
                            diffarray[j,i] = (diffarray[j+1,i-1]-diffarray[j,i-1])/(duplxdata[j+i]-duplxdata[j]);
                    end
                end


                diffarray[1,duplxlen-1] = (diffarray[2,2*xlen-2]-diffarray[1,2*xlen-2])/(duplxdata[2*xlen]-duplxdata[1]);

                coefficients = zeros(2*xlen);

                coefficients[1] = ydata[1];

                for i=2:2*xlen
                    coefficients[i]=diffarray[1,i-1];
                end

                return coefficients
end # function getcoefficients
	    
""" Interpolate data or function using Hermite Polynomials by determining the
proper coefficients.

Usage:
	interpolate(xdata,ydata,interolationpoints,dydata=dydx)

Inputs:
	xdata::Vector{Real} independent data
	ydata::Vector{Real} depedent data
	interpoints::Vector{Real} interpolation points requested
	dydata::Vector{T}=Vector{Real}(undef,1)) derivative of ydata w.r.t xdata

Outputs:
	interpoints::Vector{Real}
	hermitepoly::Array{Float64,1} the hermite polynomial

Notes: Current version requires derivative data.

Example:
     x = collect(range(0,stop=pi/2,length=30));

     y = sin.(x);

     dydx = cos.(x);

     x_interpolate = collect(range(0,stop=pi/2,length=90));

     \_,y\_interpolate = interpolate(x,y,x_interpolate,dydata=dxdy)


""" function interpolate(xdata::Vector{T},
                         ydata::Vector{T},
                         interpoints::Vector{T};
                         dydata::Vector{T}=Vector{Real}(undef,1)) where T <: Real

                         coeffs = getcoefficients(xdata,ydata,dydata=dydata);
                         xlen = length(xdata);
                         interpointslen = length(interpoints);

                         duplxdata = zeros(2*xlen);
                         for i=1:xlen
                             duplxdata[2*i-1] = xdata[i];
                             duplxdata[2*i] = xdata[i];
                         end

                         hermitepoly = zeros(Float64,interpointslen);
                         diffarray = zeros(Float64,2*xlen);
                         for i=1:interpointslen
                             diffarray[1] = 1.00e0;
                             hermitepoly[i] = coeffs[1];
                             for j=2:2*xlen
                                 diffarray[j] = (interpoints[i]-duplxdata[j-1])*diffarray[j-1];
                                 hermitepoly[i] += coeffs[j]*diffarray[j];
                             end
                         end
                         return(interpoints,hermitepoly)
end #function interpolate

end #HermiteInterpolate
