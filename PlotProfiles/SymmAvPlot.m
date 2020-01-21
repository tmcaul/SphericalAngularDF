%p=Sphere.psi_pats_experiment_NMF_2{3,2};

function SymmAvPlot(p,varargin)
    x = linspace(-1,1,1000);
    input = acos(x)./degree;

    if size(p,2) ~= 1
        for i=1:size(p,2)
            y(i,:) =(p{i}.eval(x));
        end
    y_av=mean(y,1);
    else
    y_av=p.eval(x);
    end
    
    plot(input,y_av,varargin{:})
    
end