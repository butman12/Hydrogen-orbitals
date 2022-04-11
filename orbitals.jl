using Plots
theta = LinRange(0,2*pi,501)
e=MathConstants.e
fi = LinRange(0,2*pi,501)
g(x) = x^2
r1 = 0.53 * 10^-10 #1st bohr radius
xi = yi = LinRange(0,25*r1,501) #Domain - range from the mass center
"""
functions whic name starts with "p" are polar functions, describing behaviour in respect to theta and fi angle

functions whic name starts with "r" are radial funkctions, describing behaviour in respect to distance feom the cemter of mass

Used equations are directly taken from .pdf file where topic is described in wider terms, therefore
authors see no need to explain them here
"""
p00(t,p) = (abs(1/(2*pi)))^2
p10(t,p) = (abs(((3/(4*pi))^(1/2))*cos(t)))^2
p11(t,p) = (abs(((3/(8*pi))^(1/2))*sin(t)*e^(im*p)))^2
p20(t,p) = (abs(((5/(16*pi))^(1/2))*((3*(cos(t))^2)-1)))^2
p21(t,p) = (abs(((15/(8*pi))^(1/2))*cos(t)*sin(t)*e^(im*p)))^2
p22(t,p) = (abs(((15/(32*pi))^(1/2))*sin(t)*sin(t)*e^(2im*p)))^2
r10(x) = ((1/r1))^(3/2) * 2*e^(-x/r1)
r20(x) = (1/(2*r1))^(3/2) *2*(1-(x/(2*r1)))*e^(-x/(2*r1))
r21(x) = (1/(2*r1))^(3/2) *(x/((3^(1/2))*r1))*e^(-x/(2*r1))
r30(x) = (1/(3*r1))^(3/2) * 2*(1-((2*x)/(3*r1))+(2/27)*(x/r1)^2)*e^(-x/(3*r1))
r31(x) = (1/(3*r1))^(3/2) * 4/3*(2^(1/2))/3 *(x/r1) * (1-(x/(6*r1)))*e^(-x/(3*r1))
r32(x) = (1/(3*r1))^(3/2) * (2*2^(1/2))/(27*5^(1/2)) * ((x/r1)^2)*e^(-x/(3*r1))


"""
to draw a graph call "graph function" with arguments in order below:
1) type of function used:
    -"radial"
    -"radialprob" - radial probability function
    -"polar" - polar probability function
    -"mixed" - probability function (product of polar and radial probabilities)
- quantum numer N - max 3
- quantum number L - not greater than N-1
- quantum numer M - not greater than L (since negative M describes almost identical situation as its
positive counterpart, all negative M's will be interpreted as positive.)
"""
function graph(type, n, l, m)
    if type == "polar"
        if l == 0
            data = p00.(theta,fi)
            plot(theta, data, proj=:polar, title = "l=0, m=0", label =:none)
        elseif l == 1
            if m == 0
                data = p10.(theta,fi)
                plot(theta, data, proj=:polar, title = "l=1, m=0", label =:none)
            elseif m == 1 || m == -1
                data = p11.(theta,fi)
                plot(theta, data, proj=:polar, title = "l=1, m=1", label =:none)
            end
        elseif l == 2
            if m == 0
                data = p20.(theta,fi)
                plot(theta, data, proj=:polar, title = "l=2, m=0", label =:none)
            elseif m == 1 || m == -1
                data = p21.(theta,fi)
                plot(theta, data, proj=:polar, title = "l=2, m=1", label =:none)
            elseif m == 2 || m == -2
                data = p22.(theta,fi)
                plot(theta, data, proj=:polar, title = "l=2, m=2", label =:none)
            end
        end

    elseif type == "radialprob"
        if n == 1
            data = [r10(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
            plot(xi, data, title = "n=1, l=0", label =:none)
            xlabel!("distance form the center of mass [m]")
            ylabel!("probability of elecron's position")
        elseif n == 2
            if l == 0
                data = [r20(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
                plot(xi, data, title = "n=2, l=0", label =:none)
                xlabel!("distance form the center of mass [m]")
                ylabel!("probability of elecron's position")
            elseif l == 1
                data = [r21(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
                plot(xi, data, title = "n=2, l=1", label =:none)
                xlabel!("distance form the center of mass [m]")
                ylabel!("probability of elecron's position")
            end
        elseif n == 3
            if l == 0
                data = [r30(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
                plot(xi, data, title = "n=3, l=0", label =:none)
                xlabel!("distance form the center of mass [m]")
                ylabel!("probability of elecron's position")
            elseif l == 1
                data = [r31(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
                plot(xi, data, title = "n=3, l=1", label =:none)
                xlabel!("distance form the center of mass [m]")
                ylabel!("probability of elecron's position")
            elseif l == 2
                data = [r32(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
                plot(xi, data, title = "n=3, l=2", label =:none)
                xlabel!("distance form the center of mass [m]")
                ylabel!("probability of elecron's position")
            end
        end

    elseif type == "radial"
        if n == 1
            data = r10.(xi)
            plot(xi, data, title = "n=1, l=0", label =:none)
            xlabel!("distance form the center of mass [m]")
        elseif n == 2
            if l == 0
                data = r20.(xi)
                plot(xi, data, title = "n=2, l=0", label =:none)
                xlabel!("distance form the center of mass [m]")
            elseif l == 1
                data = r21.(xi)
                plot(xi, data, title = "n=2, l=1", label =:none)
                xlabel!("distance form the center of mass [m]")
            end
        elseif n == 3
            if l == 0
                data = r30.(xi)
                plot(xi, data, title = "n=3, l=0", label =:none)
                xlabel!("distance form the center of mass [m]")
            elseif l == 1
                data = r31.(xi)
                plot(xi, data, title = "n=3, l=1", label =:none)
                xlabel!("distance form the center of mass [m]")
            elseif l == 2
                data = r32.(xi)
                plot(xi, data, title = "n=3, l=2", label =:none)
                xlabel!("distance form the center of mass [m]")
            end
        end
    elseif type == "mixed" #datar i datap to wartosci odpowiednio prawdopodobienstwa radialnego i polarnego
        if n == 1
            datar = [r10(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
            datap = p00.(theta, fi)
            bothfunctions = [polar * radial for polar in datap, radial in datar]
            Plots.heatmap(xi,theta,bothfunctions, title = "n=1, l=0, m=0")
            xlabel!("distance form the center of mass [m]")
            ylabel!("angle [radians]")

        elseif n == 2
            if l == 0
                datar = [r20(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
                datap = p00.(theta,fi)
                bothfunctions = [polar * radial for polar in datap, radial in datar]
                Plots.heatmap(xi,theta,bothfunctions, title = "n=2, l=0, m=0")
                xlabel!("distance form the center of mass [m]")
                ylabel!("angle [radians]")
            elseif l == 1
                datar = [r21(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
                if m == 0
                    datap = p10.(theta,fi)
                    bothfunctions = [polar * radial for polar in datap, radial in datar]
                    Plots.heatmap(xi,theta,bothfunctions, title = "n=2, l=1, m=0")
                    xlabel!("distance form the center of mass [m]")
                    ylabel!("angle [radians]")
                elseif m == 1 || m == -1
                    datap = p11.(theta,fi)
                    bothfunctions = [polar * radial for polar in datap, radial in datar]
                    Plots.heatmap(xi,theta,bothfunctions, title = "n=2, l=1, m=1")
                    xlabel!("distance form the center of mass [m]")
                    ylabel!("angle [radians]")
                end
            end
        elseif n == 3
                if l == 0
                datar = [r30(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
                datap = p00.(theta,fi)
                bothfunctions = [polar * radial for polar in datap, radial in datar]
                Plots.heatmap(xi,theta,bothfunctions, title = "n=3, l=0, m=0")
                xlabel!("distance form the center of mass [m]")
                ylabel!("angle [radians]")
            elseif l == 1
                datar = [r31(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
                if m == 0
                    datap = p10.(theta,fi)
                    bothfunctions = [polar * radial for polar in datap, radial in datar]
                    Plots.heatmap(xi,theta,bothfunctions, title = "n=3, l=1, m=0")
                    xlabel!("distance form the center of mass [m]")
                    ylabel!("angle [radians]")
                elseif m == 1 || m == -1
                    datap = p11.(theta,fi)
                    bothfunctions = [polar * radial for polar in datap, radial in datar]
                    Plots.heatmap(xi,theta,bothfunctions, title = "n=3, l=1, m=1")
                    xlabel!("distance form the center of mass [m]")
                    ylabel!("angle [radians]")
                end
            elseif l == 2
                datar = [r32(x)^2 * g(x)/10^10/(0.53)^-1 for x in xi]
                if m == 0
                    datap = p20.(theta,fi)
                    bothfunctions = [polar * radial for polar in datap, radial in datar]
                    Plots.heatmap(xi,theta,bothfunctions, title = "n=3, l=2, m=0")
                    xlabel!("distance form the center of mass [m]")
                    ylabel!("angle [radians]")
                elseif m == 1 || m == -1
                    datap = p21.(theta,fi)
                    bothfunctions = [polar * radial for polar in datap, radial in datar]
                    Plots.heatmap(xi,theta,bothfunctions, title = "n=3, l=2, m=1")
                    xlabel!("distance form the center of mass [m]")
                    ylabel!("angle [radians]")
                elseif m == 2 || m == -2
                    datap = p22.(theta,fi)
                    bothfunctions = [polar * radial for polar in datap, radial in datar]
                    Plots.heatmap(xi,theta,bothfunctions, title = "n=2, l=2, m=2")
                    xlabel!("distance form the center of mass [m]")
                    ylabel!("angle [radians]")
                end
            end
        end
    end
end



graph("mixed",3,2,2)
