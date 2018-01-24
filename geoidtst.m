function [geoid,msg] = geoidtst(geoid)
if nargout ~= 0;  msg = [];  end
%  Test inputs
if nargin ~= 1;  error('Incorrect number of arguments');   end
%  Ensure a real input
if ~isreal(geoid)
	  warning('Imaginary part of complex GEOID input ignored')
      geoid = real(geoid);
end
%  Geoid vector tests
if isstr(geoid)
      msg = 'Geoid vector must have 1 or 2 elements';
	  if nargout < 2;  error(msg);  end
	  return
elseif max(size(geoid)) == 1
	  geoid = [geoid 0];
elseif ~isequal(sort(size(geoid)),[1 2])
      msg = 'Geoid vector must have 1 or 2 elements';
	  if nargout < 2;  error(msg);  end
	  return
elseif (geoid(2) >= 1) | (geoid(2) < 0)
      msg = 'Geoid eccentricity must be in [0,1]';
	  if nargout < 2;  error(msg);  end
	  return
end