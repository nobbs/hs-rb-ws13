function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[1 0.99999999999999989 1]);
set(ax,'PlotBoxAspectRatio',[3.0000000000000004 2.0000000000000004 1]);
set(ax,'XLimMode','auto');
set(ax,'YLimMode','auto');
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');

% Geometry description:
pderect([0 0.33333333333333331 0.33333333333333331 0],'R1');
pderect([0.33333333333333331 0.66666666666666663 0.33333333333333331 0],'R2');
pderect([0.66666666666666663 1 0.33333333333333331 0],'R3');
pderect([0 0.33333333333333331 0.66666666666666663 0.33333333333333331],'R4');
pderect([0.33333333333333331 0.66666666666666663 0.66666666666666663 0.33333333333333331],'R5');
pderect([0.66666666666666663 1 0.66666666666666663 0.33333333333333331],'R6');
pderect([0 0.33333333333333331 1 0.66666666666666663],'R7');
pderect([0.33333333333333331 0.66666666666666663 1 0.66666666666666663],'R8');
pderect([0.66666666666666663 1 1 0.66666666666666663],'R9');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1+R2+R3+R4+R5+R6+R7+R8+R9')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'10.0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0 ';...
'0.0 ';...
'10.0';...
'1.0 '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','1000','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');
