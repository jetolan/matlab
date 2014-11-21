%function print_fig(name, fig_handle, high_res, renderer)
%
% Creates pdf in vector format and high resolution png figure 
%
%  Inputs: 
%    
%   name: filename of output figures will be name.pdf and name.png
%   fig_handle: figure handle 
%   high_res  : is true, then png will have crazy high resolution
%   renderer  : either 'manual' (default) or 'auto'. Auto will
%               sometimes not make a vector, depending on the image
%
%%%%%%%%%%%%%%%%%%%%%%%%
function print_fig(name, fig_handle, high_res)

if(~exist('fig_handle', 'var')|| isempty(fig_handle))
  fig_handle=gcf;
end

if(~exist('high_res', 'var')|| isempty(high_res))
  high_res=1;
end


if high_res
  density='500'; %resolution of png in pixels per inch
else
  density='200';
end

set(fig_handle,'PaperPositionMode','auto');  %set printing to window size
set(fig_handle, 'Renderer', 'painters')
set(fig_handle, 'RendererMode', 'manual')
print(fig_handle, '-depsc2', '-Painters',name); %need painters for vector
fix_lines([name '.eps'])
system_safe(['epstopdf ' name '.eps']);
%system_safe(['rm ' name '.eps']);  %maybe we want to keep eps for 
                                    %map images, which look better if ...
				    %you convert to pdf using an ...
				    %external program
				      

system_safe(['convert -density ' density  ' -units pixelsperinch ' name '.pdf ' name '.png']);

return
