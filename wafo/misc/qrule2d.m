function [bpx,bpy,wfxy] = qrule2d (nx,ny)
%QRULE2D Computes nodes and weights for Gaussian quadratures 
%
% CALL:  [bpx,bpy,wfxy]=qrule2d(nx,ny)
%  
%  bpx, bpy = base points (abscissas)
%  wfxy     = weight factors
%  nx,ny    = number of base points (abscissas) in each direction.
%
%  See also  gaussq2d


[bpxv,wfxv]=qrule(nx);
[bpyv,wfyv]=qrule(ny);
[bpx,bpy]=meshgrid(bpxv,bpyv);
[wfx,wfy]=meshgrid(wfxv,wfyv);
wfxy=wfx.*wfy;
