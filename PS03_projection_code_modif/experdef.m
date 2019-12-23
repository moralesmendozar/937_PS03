function allexpers = experdef()

%%% Define parameters and grids

%========================================================================================
%   set parameters
%========================================================================================

params=struct;
% AR(1) productivity shock
params.a=1;
params.alower=.9;
params.sig_a=.01;
%params.sig_a=.0001;
params.rho_a=0.5;
params.N_a=5; % Number of discrete states

% preferences
params.gamma=2;
params.betaF=0.95;
params.betaG=0.98;

% technology
params.alpha=0.7;
params.Kbar=1;
params.phi=0;

% credit
params.theta=0.9;
params.use_min_q = true;

%========================================================================================
%   grid definitions below
%========================================================================================

% grid for bench
grids.benchgrid.kF={ linspace(0.06, 0.8, 8) };
grids.benchgrid.levF={ linspace(0.75,1.0,10) };

% grid for theta=0.875
grids.theta875grid.kF = { linspace(0.01, 0.99, 10) };
grids.theta875grid.levF = { linspace(0.75,0.95,12) };

% grid for alpha=0.90
grids.alpha90grid.kF = { linspace(0.01, 0.99, 10) };
grids.alpha90grid.levF = { linspace(0.75,1.00,12) };


defaultgrid = grids.benchgrid;

%========================================================================================
%   grid definitions above
%========================================================================================

% Define experiments as modifications to base
expers = {'bench',{};
          'bench_2',{};
		  'theta875',{'theta',0.875; 'grid','theta875grid'};
		  'alpha90',{'alpha',0.9; 'alower',0.85; 'grid','alpha90grid'};
				};

%% Create each experiment definition (fully dynamic below this line)
N_EXP=size(expers,1); % = 3 in this case :)

gridvarnames = fieldnames(defaultgrid)';
N_VAR = numel(gridvarnames);

paramnames = fieldnames(params)';
paramtypes = struct2cell(structfun(@(x)class(x),params,'UniformOutput',false))';
N_PARAM = numel(paramnames);

% Write list of defined experiments to file
fid = fopen([mfilename,'.txt'],'w');
for i=1:N_EXP
   fprintf(fid,expers{i,1});
   fprintf(fid,'\n');
end
fclose(fid);

% Create experiment table
expertable = table('Size', [N_EXP,1+N_VAR+N_PARAM], ...
	'VariableTypes', [{'string'}, ...
					  repmat({'cell'},1,N_VAR), ...
					  paramtypes] ,...
	'VariableNames', [{'exper'}, gridvarnames, paramnames]);

% Fill with default values
expertable.exper(:,1) = expers(:,1);
expertable(:,2:N_VAR+1) = repmat( struct2cell( defaultgrid )', N_EXP, 1 );
expertable(:,N_VAR+2:end) = repmat( struct2cell( params )', N_EXP, 1);
expertable.Properties.RowNames = expertable.exper;

% Process modifications
for i=1:N_EXP
   changes = expers{i,2};
   for j=1:size(changes,1)
      varname=changes{j,1};
      if strcmp(varname,'grid')
          newgrid = struct2cell(grids.(changes{j,2}))';
          expertable(i,2:N_VAR+1) = newgrid;
      else
          if ischar(changes{j,2})
            expertable.(varname)(i,:)=changes(j,2);
          else
            expertable.(varname)(i,:)=changes{j,2};    
          end
      end
   end
end

% Create a struct of structs
allexpers=cell2struct(expertable.Properties.RowNames,expertable.Properties.RowNames);
for i=1:size(expertable,1)
    s.grids = table2struct(expertable(i,2:N_VAR+1));
    s.params = table2struct(expertable(i,N_VAR+2:end));
    name = expertable.Properties.RowNames{i};
    allexpers.(name) = s;
end