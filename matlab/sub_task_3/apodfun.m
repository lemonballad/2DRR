function [apod] = apodfun(x,x0,gL,gH,wid,invert)
%Generates an asymmetric curve (bandpass filter) decaying on each side with a flat top
%
%if invert=1, function will invert to be a notch filter instead of a
%bandpass filter
%
%x: A row or column vector, the domain of the function.
%x0: The desired center of the function.
%gL: The decay on the low side of the function.
%gH: The decay on the high side of the function.
%wid: The TOTAL width of the function's flat top, such that the distance
%from x0 to the portion of decay is (wid/2).
%
%
locent=x0-(wid/2);
hicent=x0+(wid/2);
apod(length(x))=0;
fu(length(x))=0;

for kk=1:length(x)
  if x(kk)<locent %REGION 1
    fu(kk)=exp(-0.5*(x(kk)-locent)^2/gL^2);
  elseif x(kk)>hicent %REGION 3
    fu(kk)=exp(-0.5*(x(kk)-hicent)^2/gH^2);
  else %REGION 2
      fu(kk)=1;
  end
end

if invert==1
    apod=1-fu;
else
    apod=fu;
end

return;