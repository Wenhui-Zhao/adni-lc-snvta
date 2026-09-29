function figure4_plsc(mode, results_root, toolbox_root, task_ids)
% Figure 4 scalar PLSC using the executed v1d McIntosh MATLAB engine.
% Inputs: .work/plsc_inputs/<LC|SNVTA>/<lme_intercept|lme_slope>/
%         X.csv, Y.csv, rows.csv, scaling.csv; X/Y start with ordered RID.
% Four direct scalar tasks use fixed study resampling settings.
if nargin < 4, task_ids = []; end
assert(nargin >= 3, 'Figure4:Arguments', 'Mode, results root and toolbox root required.');
assert(~isempty(results_root),'Figure4:Output','Explicit output root required.');
out = fullfile(results_root, '.work', 'plsc');
switch char(mode)
    case 'prepare'
        prepare_inputs(results_root, toolbox_root, false);
    case 'task'
        assert(isscalar(task_ids), 'Figure4:TaskID', 'Exactly one task ID required.');
        run_task(out, toolbox_root, task_ids);
    case 'finalize'
        if isempty(task_ids), task_ids=selected_ids(task_table(false)); end
        finalize_outputs(out, task_ids);
    otherwise
        error('Figure4:Mode', 'Unknown PLSC mode: %s', mode);
end
end

function T = task_table(synthetic)
ids = (1:4)';
branch = ["LC";"LC";"SNVTA";"SNVTA"];
model = repmat(["lme_intercept";"lme_slope"],2,1); nucleus=branch;
variant=repmat("direct",4,1);seed=repmat([20260721;20260722],2,1);
nperm=repmat(5000,4,1);nboot=nperm;profile=repmat("study",4,1);
comparison=repmat("current_main_replication",4,1);
T=table(ids,branch,model,nucleus,variant,seed,nperm,nboot,profile,comparison, ...
'VariableNames',{'task_id','branch','model','nucleus','variant','seed','num_permutations','num_bootstraps','resampling_profile','comparison_class'});
end

function prepare_inputs(root, toolbox_root, synthetic)
out = fullfile(root,'.work','plsc');
if exist(fullfile(out,'prepared_manifest.mat'),'file')
    P=validate_prepared(out,toolbox_root);
    assert(P.manifest.synthetic==synthetic,'Figure4:Profile','Synthetic/study preparation mismatch.');
    refresh_pending(out,P.tasks);
    fprintf('Verified unchanged prepared scalar PLSC bundle; pending IDs refreshed.\n');
    return;
end
if ~exist(out,'dir'), mkdir(out); end
mkdir(fullfile(out,'inputs')); mkdir(fullfile(out,'tasks'));
tasks = task_table(synthetic);
toolbox = toolbox_identity(toolbox_root);
source_hash = sha256_file(mfilename('fullpath') + ".m");
models = {'lme_intercept','lme_slope'}; nuclei = {'LC','SNVTA'};
input_records = struct('path',{},'sha256',{});
source_records = struct('path',{},'sha256',{});
regions_path=fullfile(root,'.work','plsc_inputs','regions.csv');
if ~synthetic
    assert(exist(regions_path,'file')==2,'Figure4:Regions','Authoritative primary145 region table missing.');
    regions=readtable(regions_path,'FileType','text','Delimiter',',');
    assert(height(regions)==145 && ismember('variable_column',regions.Properties.VariableNames) && ...
        ~ismember('feature',regions.Properties.VariableNames), ...
        'Figure4:Regions','Region table requires the authoritative primary145 variable_column schema.');
    source_records(end+1)=struct('path',regions_path,'sha256',sha256_file(regions_path));
end
for p=1:2
    B = cell(1,2);
    for n=1:2
        base = fullfile(root,'.work','plsc_inputs',nuclei{n},models{p});
        B{n} = read_direct(base, synthetic);
        source_records=[source_records,B{n}.source_files]; %#ok<AGROW>
        if ~synthetic
            assert(isequal(string(B{n}.features(:)),string(regions.variable_column)), ...
                'Figure4:FeatureOrder','Regional matrix order differs from authoritative primary145.');
            B{n}.region_definition_path=regions_path;
            B{n}.region_definition_sha256=sha256_file(regions_path);
        end
        bundle = B{n}; %#ok<NASGU>
        fp = fullfile(out,'inputs',sprintf('%s_%s.mat',nuclei{n},models{p}));
        save(fp,'bundle','-v7.3');
        input_records(end+1) = struct('path',fp,'sha256',sha256_file(fp)); %#ok<AGROW>
    end

end
manifest = struct('schema','figure4_scalar_plsc_v1','execution_scope','main','synthetic',synthetic, ...
    'source_sha256',source_hash,'toolbox',toolbox,'inputs',input_records,'trajectory_sources',source_records, ...
    'runtime_version',version,'runtime_release',version('-release'), ...
    'sign_rule','largest_absolute_behavior_salience_positive', ...
    'seed_rule','parameter_specific_original_seed', ...
    'model_signature',getenv('MODEL_SIGNATURE'),'model_profile',runtime_profile(), ...
    'bootstrap_rng_policy',bootstrap_rng_identity(toolbox_root), ...
    'scientific_scope','main');
writetable(tasks, fullfile(out,'tasks.csv'),'FileType','text','Delimiter',',');
manifest.tasks_sha256 = sha256_file(fullfile(out,'tasks.csv'));
save(fullfile(out,'prepared_manifest.mat'),'manifest','tasks','-v7');
write_json(fullfile(out,'prepared_manifest.json'),manifest);
refresh_pending(out,tasks);
fprintf('Prepared %d scalar PLSC tasks in %d input bundles; profile=%s; scope=%s.\n',height(tasks),numel(input_records),tasks.resampling_profile(1),'main');
end

function B = read_direct(base, synthetic)
names = {'X.csv','Y.csv','rows.csv','scaling.csv'};
provenance = struct('path',{},'sha256',{});
for k=1:numel(names)
    fp=fullfile(base,names{k});
    assert(exist(fp,'file')==2,'Figure4:MissingInput','Required trajectory input missing: %s',fp);
    provenance(k)=struct('path',fp,'sha256',sha256_file(fp));
end
xh=csv_header(provenance(1).path);yh=csv_header(provenance(2).path);rh=csv_header(provenance(3).path);
assert(numel(unique(xh))==numel(xh),'Figure4:FeatureID','Repeated regional column names.');
XT=readtable(provenance(1).path,'FileType','text','Delimiter',',');
YT=readtable(provenance(2).path,'FileType','text','Delimiter',',');
rows=readtable(provenance(3).path,'FileType','text','Delimiter',',');
assert(isequal(xh,XT.Properties.VariableNames),'Figure4:FeatureID','Regional header was altered during import.');
assert(isequal(yh,YT.Properties.VariableNames) && isequal(rh,rows.Properties.VariableNames), ...
    'Figure4:Rows','Behavior or row header was altered during import.');
assert(strcmp(XT.Properties.VariableNames{1},'RID') && strcmp(YT.Properties.VariableNames{1},'RID'), ...
    'Figure4:Rows','X/Y must begin with RID.');
assert(ismember('RID',rows.Properties.VariableNames),'Figure4:Rows','rows.csv lacks RID.');
XT.RID=ordered_numeric_ids(XT.RID,height(XT));
YT.RID=ordered_numeric_ids(YT.RID,height(YT));
rows.RID=ordered_numeric_ids(rows.RID,height(rows));
assert(isequal(XT.RID,YT.RID) && isequal(XT.RID,rows.RID),'Figure4:Rows','Ordered X/Y/rows identities differ.');
assert(width(XT)==146 && width(YT)==2,'Figure4:Dimension','Exactly 145 regional columns and scalar Y required.');
X=declared_numeric_columns(XT(:,2:end)); Y=declared_numeric_columns(YT(:,2:end));
features=XT.Properties.VariableNames(2:end);
tokens=regexp(features,'^h_muse_volume_([0-9]+)$','tokens','once');
assert(all(~cellfun(@isempty,tokens)),'Figure4:FeatureID','Unexpected primary145 feature spelling.');
region_id=cellfun(@(x)str2double(x{1}),tokens)';
assert(numel(unique(region_id))==145,'Figure4:FeatureID','Repeated atlas region IDs.');
assert(all(isfinite(X(:))) && all(isfinite(Y(:))),'Figure4:Nonfinite','Trajectory final matrices must be fully finite.');
assert(size(X,1)>=10,'Figure4:Sample','PLSC requires at least ten eligible participants.');
assert(all(std(X,0,1)>0) && std(Y)>0,'Figure4:Scale','Constant regional or behavior parameter.');
B=struct('X',X,'Y',Y,'RID',XT.RID,'rows',rows,'features',{features}, ...
    'region_id',region_id,'behavior',YT.Properties.VariableNames{2}, ...
    'source_files',provenance,'trajectory_scaling_file',provenance(4).path, ...
    'reference_population','complete_finite_X145_and_bilateral_Y_participants_in_nucleus_parameter_branch', ...
    'synthetic',synthetic);
end

function ids = ordered_numeric_ids(value,n)
% R's precision-preserving CSV writer quotes character-encoded numeric IDs.
% Accept exact integer IDs only; never evaluate text or repair row order.
if isnumeric(value)
    assert(isreal(value) && iscolumn(value) && numel(value)==n, ...
        'Figure4:Rows','RID must be a real ordered column.');
    if isinteger(value)
        assert(all(value<=cast(flintmax,'like',value)) && ...
            all(value>=cast(-flintmax,'like',value)), ...
            'Figure4:Rows','RID exceeds exact double integer representation.');
    elseif isa(value,'single')
        assert(all(abs(value)<=flintmax('single')), ...
            'Figure4:Rows','RID exceeds exact single integer representation.');
    end
    ids=double(value);
elseif isstring(value) || iscellstr(value)
    value=string(value);
    assert(iscolumn(value) && numel(value)==n && ~any(ismissing(value)), ...
        'Figure4:Rows','RID text must be a complete ordered column.');
    text_ids=cellstr(value);
    for k=1:n
        s=text_ids{k};
        assert(~isempty(regexp(s,'^[+-]?[0-9]+(?:\.0+)?$','once')), ...
            'Figure4:Rows','RID requires decimal integer text, without blanks or prefixes.');
        digits=regexprep(regexprep(s,'^[+-]',''),'\.0+$','');
        digits=regexprep(digits,'^0+','');
        if isempty(digits),digits='0';end
        limit='9007199254740992';
        assert(numel(digits)<numel(limit) || ...
            (numel(digits)==numel(limit) && isequal(sort({digits,limit}),{digits,limit})), ...
            'Figure4:Rows','RID text exceeds exact double integer representation.');
    end
    ids=str2double(value);
else
    error('Figure4:Rows','Unsupported RID representation.');
end
assert(all(isfinite(ids)) && all(abs(ids)<=flintmax) && all(ids==fix(ids)) && ...
    numel(unique(ids))==n,'Figure4:Rows', ...
    'RID must contain distinct finite exactly represented integers; conversion collisions are forbidden.');
end

function names = csv_header(fp)
% Producer names are simple identifiers; readtable must not silently repair them.
fid=fopen(fp,'r');assert(fid>=0,'Figure4:Unreadable','Cannot read matrix header.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
line=fgetl(fid);assert(ischar(line),'Figure4:Header','Empty matrix file.');
names=strsplit(line,',','CollapseDelimiters',false);
for k=1:numel(names)
    s=names{k};
    if numel(s)>=2 && s(1)=='"' && s(end)=='"',s=s(2:end-1);end
    assert(~isempty(s) && ~any(s=='"'),'Figure4:Header','Malformed matrix header.');
    names{k}=s;
end
end

function A = declared_numeric_columns(T)
% atomic_csv preserves double precision as quoted %.17g text. Convert only
% the declared measurement fields; retain original rows, columns and values.
A=zeros(height(T),width(T));
for j=1:width(T)
    value=T.(T.Properties.VariableNames{j});
    if isnumeric(value)
        assert(isreal(value) && iscolumn(value) && numel(value)==height(T), ...
            'Figure4:NumericInput','Measurement must be a real numeric column.');
        v=double(value);
        if isinteger(value)
            assert(isequal(cast(v,'like',value),value), ...
                'Figure4:NumericInput','Integer measurement loses precision in double conversion.');
        end
    elseif isstring(value) || iscellstr(value)
        value=string(value);
        assert(iscolumn(value) && numel(value)==height(T) && ~any(ismissing(value)), ...
            'Figure4:NumericInput','Numeric measurement text must be complete.');
        tokens=cellstr(value);
        % This is the producer's decimal/scientific syntax, not MATLAB code.
        valid=cellfun(@(x)~isempty(regexp(x,'^-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:e[+-][0-9]+)?$','once')),tokens);
        assert(all(valid),'Figure4:NumericInput','Malformed numeric measurement token.');
        v=str2double(value);
        mantissa=regexprep(tokens,'e[+-][0-9]+$','');
        nonzero=cellfun(@(x)~isempty(regexp(x,'[1-9]','once')),mantissa);
        assert(~any(v==0 & nonzero),'Figure4:NumericInput','Nonzero measurement underflows to zero.');
    else
        error('Figure4:NumericInput','Unsupported numeric measurement representation.');
    end
    assert(all(isfinite(v)) && isreal(v), ...
        'Figure4:NumericInput','Measurements must remain finite and real after conversion.');
    A(:,j)=v;
end
end

function run_task(out, toolbox_root, task_id)
P=validate_prepared(out,toolbox_root);
T=P.tasks(P.tasks.task_id==task_id,:);
assert(height(T)==1,'Figure4:TaskID','Unsupported scientific task ID.');
assert(maxNumCompThreads==1,'Figure4:Threads','One MATLAB computation thread is required.');
dest=fullfile(out,'tasks',sprintf('%04d',task_id));
if ~exist(dest,'dir'), mkdir(dest); end
work=out;
if exist(fullfile(work,sprintf('task_%04d.json',task_id)),'file')
    validate_task(out,task_id);
    fprintf('Verified completed PLSC task %d; no repeat fit or export.\n',task_id);
    return;
end
[B,X,Y,input_path,behavior]=task_input(out,T);
x_center=mean(X,1); y_center=mean(Y,1);
X=bsxfun(@minus,X,x_center); Y=bsxfun(@minus,Y,y_center);
if exist(fullfile(dest,'result.mat'),'file')
    saved=load(fullfile(dest,'result.mat'),'result','receipt','transformation');
    assert(strcmp(saved.receipt.source_sha256,P.manifest.source_sha256) && ...
        strcmp(saved.receipt.model_signature,P.manifest.model_signature) && ...
        strcmp(saved.receipt.input_sha256,sha256_file(input_path)), ...
        'Figure4:ResumeIdentity','Saved fit identity differs; no replacement permitted.');
    % No success receipt exists: recover exports from the saved fit, never refit.
    export_task(dest,B,X,Y,saved.result,plsc_bootstrap_ratios(saved.result,1),T,saved.receipt,saved.transformation);
    write_json(fullfile(dest,'runtime_receipt.json'),saved.receipt);
    write_task_receipt(out,task_id,saved.receipt);
    fprintf('Verified completed PLSC task %d; no repeat fit.\n',task_id);
    return;
end
[path_cleanup, resolved_runtime]=scoped_engine(toolbox_root); %#ok<ASGLU>

opt=struct('method',3,'num_perm',T.num_permutations,'num_boot',T.num_bootstraps, ...
    'clim',95,'stacked_behavdata',Y,'cormode',0,'boot_type','strat');
lastwarn(''); rng(double(T.seed),'twister');
rng_before=rng; start=tic;
result=pls_analysis({X},size(X,1),1,opt);
elapsed_seconds=toc(start); rng_after=rng; [warning_message,warning_id]=lastwarn;
bsr=plsc_bootstrap_ratios(result,1);
resampling=resampling_diagnostics(result,T.num_permutations,T.num_bootstraps);
assert(numel(result.s)==1 && size(result.u,2)==1 && size(result.v,1)==1, ...
    'Figure4:ScalarLV','Scalar PLSC must return one latent variable.');
assert(isfield(result,'usc') && isfield(result,'vsc'),'Figure4:Scores','Original toolbox did not return participant scores.');
% The original engine casts input to single on MATLAB >=7; reproduce that operation.
native_score=cast(X,'like',result.u)*result.u(:,1);
score_error=max(abs(double(result.usc(:,1))-double(native_score)));
assert(score_error==0,'Figure4:Scores','Native all145 score does not reproduce the toolbox score.');
assert(isfield(result,'perm_result') && isfield(result.perm_result,'sprob'), ...
    'Figure4:Inference','Permutation result missing.');
assert(isfinite(result.perm_result.sprob(1)) && result.perm_result.sprob(1)>=0 && result.perm_result.sprob(1)<=1, ...
    'Figure4:Inference','Invalid global permutation probability.');
orientation_factor=1; if double(result.v(1))<0, orientation_factor=-1; end
receipt=struct('task_id',task_id,'branch',char(T.branch),'model',char(T.model), ...
    'model_signature',P.manifest.model_signature,'model_profile',P.manifest.model_profile, ...
    'resampling_profile',char(T.resampling_profile),'seed',T.seed, ...
    'num_permutations',T.num_permutations,'num_bootstraps',T.num_bootstraps, ...
    'source_sha256',P.manifest.source_sha256,'task_sha256',P.manifest.tasks_sha256, ...
    'toolbox_fingerprint',P.manifest.toolbox.fingerprint,'engine_precision',class(result.u), ...
    'cpu_threads',maxNumCompThreads,'bootstrap_rng_policy',P.manifest.bootstrap_rng_policy, ...
    'resolved_runtime',resolved_runtime,'resampling',resampling, ...
    'num_groups',1,'num_conditions',1,'group_structure','one_group_one_condition_no_site_or_diagnosis_strata', ...
    'runtime_version',version,'runtime_release',version('-release'), ...
    'host',getenv('HOSTNAME'),'job_id',getenv('JOB_ID'),'sge_task_id',getenv('SGE_TASK_ID'), ...
    'n_subjects',size(X,1),'n_features',size(X,2),'n_behaviors',size(Y,2), ...
    'input_file',input_path,'input_sha256',sha256_file(input_path), ...
    'orientation_factor',orientation_factor,'score_max_abs_error',score_error, ...
    'elapsed_seconds',elapsed_seconds,'warning_message',warning_message,'warning_id',warning_id, ...
    'status','estimated','native_estimates_replaced',false);
transformation=struct('input_file',input_path,'input_sha256',receipt.input_sha256, ...
    'x_center_before_toolbox',x_center,'y_center_before_toolbox',y_center, ...
    'feature_names',{B.features},'region_id',B.region_id,'behavior',behavior, ...
    'orientation_factor',orientation_factor,'parameter_type',char(T.model), ...
    'engine_precision',class(result.u), ...
    'brain_score_definition','native_single(prepared_X - x_center_before_toolbox) * all145(native_result.u); then apply orientation_factor', ...
    'reference_population',B.reference_population,'input_bundle_contains_prior_transforms',true);
% Native result has the original draws/order; do not store duplicate X/RID copies.
save(fullfile(dest,'result.mat'),'result','rng_before','rng_after','receipt','transformation','T','-v7.3');
export_task(dest,B,X,Y,result,bsr,T,receipt,transformation);
write_json(fullfile(dest,'runtime_receipt.json'),receipt);
write_task_receipt(out,task_id,receipt);
fprintf('Completed scalar PLSC task %d: N=%d, features=%d, status=estimated, profile=%s.\n', ...
    task_id,size(X,1),size(X,2),T.resampling_profile);
end

function P=validate_prepared(out,toolbox_root)
P=load(fullfile(out,'prepared_manifest.mat'),'manifest','tasks');
assert(isfield(P.manifest,'execution_scope') && strcmp(P.manifest.execution_scope,'main'), ...
    'Figure4:Scope','Prepared execution scope changed.');
assert(strcmp(P.manifest.source_sha256,sha256_file(mfilename('fullpath')+".m")), ...
    'Figure4:CodeIdentity','PLSC source differs from the prepared signed source.');
assert(strcmp(P.manifest.tasks_sha256,sha256_file(fullfile(out,'tasks.csv'))), ...
    'Figure4:TaskIdentity','PLSC task specification changed.');
assert(isequal(P.manifest.toolbox,toolbox_identity(toolbox_root)), ...
    'Figure4:ToolboxIdentity','Original toolbox identity changed.');
assert(strcmp(P.manifest.runtime_version,version),'Figure4:RuntimeIdentity','MATLAB runtime changed.');
assert(isequal(P.manifest.bootstrap_rng_policy,bootstrap_rng_identity(toolbox_root)), ...
    'Figure4:BootstrapRNGIdentity','Bootstrap RNG policy or helper identity changed.');
assert(strcmp(P.manifest.model_signature,getenv('MODEL_SIGNATURE')) && ...
    strcmp(P.manifest.model_profile,runtime_profile()),'Figure4:Signature','Runtime profile or run signature changed.');
identities=[P.manifest.inputs,P.manifest.trajectory_sources];
for k=1:numel(identities)
    I=identities(k);
    assert(strcmp(I.sha256,sha256_file(I.path)),'Figure4:InputIdentity','Prepared PLSC input changed.');
end
end

function ids=selected_ids(T)
value=strtrim(getenv('FIGURE4_TASK_IDS'));
if isempty(value), ids=T.task_id'; return; end
tokens=regexp(value,'[, ]+','split'); ids=str2double(tokens);
assert(all(isfinite(ids)) && all(ids==fix(ids)) && numel(unique(ids))==numel(ids) && ...
    all(ismember(ids,T.task_id)),'Figure4:Selection','Unsupported or repeated PLSC task selection.');
end

function value=runtime_profile()
value=getenv('MODEL_PROFILE'); if isempty(value),value=getenv('PROFILE');end
end

function refresh_pending(out,T)
work=out; if ~exist(work,'dir'),mkdir(work);end
selected=selected_ids(T); pending=[];
for id=selected
    fp=fullfile(work,sprintf('task_%04d.json',id));
    if ~exist(fp,'file'),pending(end+1)=id;continue;end %#ok<AGROW>
    validate_task(out,id);
end
fid=fopen(fullfile(work,'n_tasks.txt'),'w');fprintf(fid,'%d\n',max(T.task_id));fclose(fid);
fid=fopen(fullfile(work,'pending_task_ids.txt'),'w');fprintf(fid,'%d\n',pending);fclose(fid);
write_json(fullfile(work,'prepared.json'),struct('model_signature',getenv('MODEL_SIGNATURE'), ...
    'source_sha256',sha256_file(mfilename('fullpath')+".m"),'n_tasks',max(T.task_id), ...
    'supported_task_ids',T.task_id','selected_task_ids',selected,'pending_task_ids',pending, ...
    'status','prepared','prepared_manifest_file',fullfile(out,'prepared_manifest.mat')));
end

function write_task_receipt(out,task_id,receipt)
work=out; if ~exist(work,'dir'),mkdir(work);end
fp=fullfile(work,sprintf('task_%04d.json',task_id));
if exist(fp,'file'),validate_task(out,task_id);return;end
receipt.result_file=fullfile(out,'tasks',sprintf('%04d',task_id),'result.mat');
receipt.result_sha256=sha256_file(receipt.result_file);
files=export_files(out,task_id);
receipt.exports=struct('path',{},'sha256',{});
for k=1:numel(files)
    receipt.exports(k)=struct('path',files{k},'sha256',sha256_file(files{k}));
end
receipt.required_files=[{receipt.result_file},files];
write_json(fp,receipt);
end

function files=export_files(out,task_id)
dest=fullfile(out,'tasks',sprintf('%04d',task_id));
names={'global.csv','regional.csv','scores.csv','behavior.csv', ...
    'runtime_receipt.json','score_transformation.json','transforms.csv'};
files=cellfun(@(x)fullfile(dest,x),names,'UniformOutput',false);
end

function R=validate_task(out,task_id)
work=out;
fp=fullfile(work,sprintf('task_%04d.json',task_id));
assert(exist(fp,'file')==2,'Figure4:Incomplete','Completed task receipt is missing.');
R=jsondecode(fileread(fp));
expected_result=fullfile(out,'tasks',sprintf('%04d',task_id),'result.mat');
assert(R.task_id==task_id && strcmp(R.status,'estimated') && ...
    strcmp(R.model_signature,getenv('MODEL_SIGNATURE')) && ...
    strcmp(R.source_sha256,sha256_file(mfilename('fullpath')+".m")) && ...
    strcmp(R.task_sha256,sha256_file(fullfile(out,'tasks.csv'))) && ...
    strcmp(R.result_file,expected_result) && strcmp(R.result_sha256,sha256_file(expected_result)), ...
    'Figure4:ResumeIdentity','Completed task identity or native result changed.');
expected=export_files(out,task_id);
assert(isfield(R,'exports') && numel(R.exports)==numel(expected), ...
    'Figure4:ExportIdentity','Completed task export catalog is incomplete.');
for k=1:numel(expected)
    assert(strcmp(R.exports(k).path,expected{k}) && ...
        strcmp(R.exports(k).sha256,sha256_file(expected{k})), ...
        'Figure4:ExportIdentity','Completed PLSC export was modified.');
end
end

function [B,X,Y,fp,behavior] = task_input(out,T)
fp=fullfile(out,'inputs',sprintf('%s_%s.mat',T.nucleus,T.model));
S=load(fp,'bundle');B=S.bundle;X=B.X;Y=B.Y;behavior=B.behavior;
end

function export_task(dest,B,X,Y,result,bsr,T,receipt,transform)
u=double(result.u(:,1)); s=double(result.s(1)); v=double(result.v(1));
factor=receipt.orientation_factor; n=size(X,1); p=size(X,2);
map_id="primary_dxbl_all__"+T.branch+"__"+T.model+"__LV1__salience";
global_result=table(T.task_id,T.branch,T.model,"LV1",map_id, ...
    double(result.perm_result.sprob(1)),double(result.perm_result.sprob(1))<.05, ...
    s,1,double(result.lvcorrs(1)),n,p,1,T.seed,T.num_permutations,T.num_bootstraps, ...
    T.resampling_profile,T.comparison_class,factor,"estimated", ...
    'VariableNames',{'task_id','branch','model','lv','map_id','permutation_p', ...
    'retained_permutation_p_lt_0_05','singular_value','fraction_crossblock_covariance', ...
    'brain_behavior_r','n_subjects','n_features','n_behaviors','seed', ...
    'num_permutations','num_bootstraps','resampling_profile','comparison_class','orientation_factor','status'});
regional=table(repmat(T.task_id,p,1),repmat(T.branch,p,1),repmat(T.model,p,1), ...
    repmat("LV1",p,1),repmat(map_id,p,1),(1:p)',B.region_id,string(B.features(:)), ...
    u,u*factor,u*s,u*s*factor,double(result.boot_result.u_se(:,1)),bsr,bsr*factor, ...
    abs(bsr)>=2.56,repmat(factor,p,1),corr(X,double(result.usc(:,1))*factor), ...
    corr(X,Y),repmat(global_result.permutation_p,p,1), ...
    'VariableNames',{'task_id','branch','model','lv','map_id','region_index','region_id','feature', ...
    'salience_raw','salience_oriented','scaled_salience_raw','scaled_salience_oriented', ...
    'bootstrap_se_scaled','bootstrap_ratio_raw','bootstrap_ratio_oriented','stable_bsr256', ...
    'orientation_factor','region_brain_score_loading','region_behavior_correlation','permutation_p'});
scores=table(repmat(T.task_id,n,1),repmat(T.branch,n,1),repmat(T.model,n,1),repmat("LV1",n,1), ...
    repmat(map_id,n,1),(1:n)',B.RID,double(result.usc(:,1)),double(result.usc(:,1))*factor, ...
    double(result.vsc(:,1)),double(result.vsc(:,1))*factor,Y,repmat(factor,n,1), ...
    'VariableNames',{'task_id','branch','model','lv','map_id','row_index','RID', ...
    'brain_score_raw','brain_score_oriented','behavior_score_raw','behavior_score_oriented', ...
    'centered_scalar_behavior','orientation_factor'});
global_result.score_evidence="descriptive_in_sample_not_independent_prediction";
behavior=table(T.task_id,T.branch,T.model,"LV1",string(transform.behavior),v,v*factor, ...
    corr(Y,double(result.usc(:,1))*factor),'VariableNames', ...
    {'task_id','branch','model','lv','behavior','salience_raw','salience_oriented','brain_score_loading'});
writetable(global_result,fullfile(dest,'global.csv'),'FileType','text','Delimiter',',');
writetable(regional,fullfile(dest,'regional.csv'),'FileType','text','Delimiter',',');
writetable(scores,fullfile(dest,'scores.csv'),'FileType','text','Delimiter',',');
writetable(behavior,fullfile(dest,'behavior.csv'),'FileType','text','Delimiter',',');
write_json(fullfile(dest,'score_transformation.json'),transform);
transform_table=table(repmat(T.branch,p,1),repmat(T.model,p,1),string(B.features(:)), ...
 transform.x_center_before_toolbox(:),repmat(transform.y_center_before_toolbox,p,1),repmat(factor,p,1), ...
 repmat(string(transform.engine_precision),p,1), ...
 'VariableNames',{'nucleus','parameter','feature','x_center_before_toolbox','y_center_before_toolbox','orientation_factor','engine_precision'});
writetable(transform_table,fullfile(dest,'transforms.csv'));
end

function finalize_outputs(out, task_ids)
S=load(fullfile(out,'prepared_manifest.mat'),'manifest');
S=validate_prepared(out,S.manifest.toolbox.root);
T=S.tasks;
if ~isempty(task_ids)
    assert(all(ismember(task_ids,T.task_id)) && numel(unique(task_ids))==numel(task_ids), ...
        'Figure4:TaskID','Invalid finalization task selection.');
    T=T(ismember(T.task_id,task_ids),:);
end
for j=1:height(T),validate_task(out,T.task_id(j));end
tables={'global','regional','scores','behavior','transforms'};
for k=1:numel(tables)
    blocks=cell(height(T),1);
    for j=1:height(T)
        dest=fullfile(out,'tasks',sprintf('%04d',T.task_id(j)));
        blocks{j}=readtable(fullfile(dest,[tables{k} '.csv']),'FileType','text','Delimiter',',');
    end
    names={'plsc_global_results.csv','plsc_regional_results.csv','plsc_participant_scores.csv', ...
        'plsc_behavior_results.csv','plsc_transforms.csv'};
    d=vertcat(blocks{:});
    if strcmp(tables{k},'global')
        d.Properties.VariableNames{strcmp(d.Properties.VariableNames,'permutation_p')}='P';
        diagnostics=d(:,{'task_id','branch','model','status'});
        for field={'warning_message','warning_id','generator','bootstrap_rng_policy','engine_precision'}
            diagnostics.(field{1})=strings(height(T),1);
        end
        for field={'actual_permutations','actual_bootstraps','validity_retry_count','effective_bootstrap_changed_columns','low_variability_bootstraps','seed','elapsed_seconds'}
            diagnostics.(field{1})=nan(height(T),1);
        end
        for row=1:height(T)
            r=jsondecode(fileread(fullfile(out,'tasks',sprintf('%04d',T.task_id(row)),'runtime_receipt.json')));
            diagnostics.warning_message(row)=string(r.warning_message);
            diagnostics.warning_id(row)=string(r.warning_id);
            diagnostics.actual_permutations(row)=r.resampling.actual_num_permutations;
            diagnostics.actual_bootstraps(row)=r.resampling.actual_num_bootstraps;
            diagnostics.validity_retry_count(row)=r.resampling.validity_retry_count;
            diagnostics.effective_bootstrap_changed_columns(row)=r.resampling.initial_to_effective_changed_columns;
            diagnostics.low_variability_bootstraps(row)=r.resampling.low_variability_bootstraps;
            diagnostics.generator(row)="twister";
            diagnostics.seed(row)=r.seed;
            diagnostics.bootstrap_rng_policy(row)=string(r.bootstrap_rng_policy.policy);
            diagnostics.engine_precision(row)=string(r.engine_precision);
            diagnostics.elapsed_seconds(row)=r.elapsed_seconds;
        end
        writetable(diagnostics,fullfile(out,'diagnostics.csv'));
        d.status=[];
    end
    root=fileparts(fileparts(out));
    if strcmp(tables{k},'scores'),dest=fullfile(root,'private','figure4');else,dest=fullfile(root,'figure4');end
    if ~exist(dest,'dir'),mkdir(dest);end
    writetable(d,fullfile(dest,names{k}));
end
expected=task_table(S.manifest.synthetic);
summary=struct('expected_task_ids',T.task_id','completed_task_count',height(T), ...
    'full_approved_grid',isequal(T.task_id,expected.task_id), ...
    'execution_scope','main','synthetic',S.manifest.synthetic, ...
    'status','completed_requested_tasks','fdr_recomputed',false, ...
    'score_use','descriptive_in_sample_not_independent_prediction', ...
    'figure5_projection','use_native_unthresholded_weights_and_saved_input_and_score_transformations', ...
    'known_incorrect_historical_bsr','compare_u_divided_by_SE_again_is_rejected_reference');
write_json(fullfile(out,'completion.json'),summary);
fprintf('Finalized %d requested scalar PLSC tasks; complete_scope=%d; synthetic=%d.\n',height(T),summary.full_approved_grid,S.manifest.synthetic);
end

function I=bootstrap_rng_identity(toolbox_root)
original=fullfile(toolbox_root,'plscmd','rri_boot_check.m');
local=fullfile(fileparts(mfilename('fullpath')),'matlab','rri_boot_check.m');
original_hash='82a8e2bf2ee427241480d5d4a2c23a262dc2673243a8b226e862a50997b7ff41';
assert(strcmp(sha256_file(original),original_hash),'Figure4:BootstrapOriginalHash', ...
    'Original rri_boot_check hash differs from the approved source.');
old=fileread(original); reset='      rng_default; rng_shuffle;';
replacement='      % rng_default; rng_shuffle; % Approved: continue the task-seeded post-permutation stream.';
assert(numel(strfind(old,reset))==1 && strcmp(fileread(local),strrep(old,reset,replacement)), ...
    'Figure4:BootstrapMinimalDiff','Local helper must suppress exactly the approved modern reset statement.');
I=struct('policy','continue_task_seeded_twister_after_original_permutations_v1', ...
    'original_helper_path',original,'original_helper_sha256',original_hash, ...
    'local_helper_path',local,'local_helper_sha256',sha256_file(local), ...
    'change','suppress_both_rng_default_and_rng_shuffle_in_modern_runtime_branch_only', ...
    'bootstrap_realization_comparison','approved_change_not_exact_historical_BSR_replication', ...
    'group_structure','one_group_one_condition_original_strat_boot_type');
end

function [cleanup, resolved]=scoped_engine(toolbox_root)
policy=bootstrap_rng_identity(toolbox_root);
prior=path; cleanup=onCleanup(@()restore_engine_path(prior));
addpath(genpath(fullfile(toolbox_root,'plscmd'))); addpath(genpath(fullfile(toolbox_root,'plsgui')));
addpath(fileparts(policy.local_helper_path),'-begin'); clear rri_boot_check;
engine=which('pls_analysis'); helper=which('rri_boot_check');
assert(strcmp(engine,fullfile(toolbox_root,'plscmd','pls_analysis.m')) && ...
    strcmp(sha256_file(engine),sha256_file(fullfile(toolbox_root,'plscmd','pls_analysis.m'))), ...
    'Figure4:EngineIdentity','Actual resolved engine differs from original pls_analysis.');
assert(strcmp(helper,policy.local_helper_path) && strcmp(sha256_file(helper),policy.local_helper_sha256), ...
    'Figure4:BootstrapResolution','Actual task helper is not the approved local copy.');
assert(get_matlab_version>=7013,'Figure4:BootstrapRuntime','Approved reset suppression requires the modern runtime branch.');
resolved=struct('engine_path',engine,'engine_sha256',sha256_file(engine), ...
    'helper_path',helper,'helper_sha256',sha256_file(helper), ...
    'matlabroot',matlabroot,'matlab_version',version,'prior_path_restored_on_success_or_error',true);
end

function restore_engine_path(prior)
path(prior); clear rri_boot_check;
end

function D=resampling_diagnostics(result,nperm,nboot)
P=result.perm_result; B=result.boot_result;
assert(P.num_perm==nperm && B.num_boot==nboot && ...
    size(P.permsamp,2)==nperm && size(B.bootsamp,2)==nboot && size(B.bootsamp_4beh,2)==nboot, ...
    'Figure4:ResamplingCount','Actual permutation/bootstrap counts differ from the declared task.');
D=struct('actual_num_permutations',double(P.num_perm),'actual_num_bootstraps',double(B.num_boot), ...
    'initial_bootstrap_columns',size(B.bootsamp,2),'effective_bootstrap_columns',size(B.bootsamp_4beh,2), ...
    'validity_retry_count',double(B.countnewtotal), ...
    'low_variability_bootstraps',double(B.num_LowVariability_behav_boots), ...
    'initial_to_effective_changed_columns',sum(any(B.bootsamp~=B.bootsamp_4beh,1)), ...
    'effective_order_field','result.boot_result.bootsamp_4beh','validity_diagnostics_field','result.boot_result.badbeh', ...
    'seed_reset_during_bootstrap',false,'native_engine_precision',class(result.u));
end

function I = toolbox_identity(root)
assert(exist(fullfile(root,'plscmd','pls_analysis.m'),'file')==2, ...
    'Figure4:Toolbox','Original McIntosh toolbox unavailable.');
files=[dir(fullfile(root,'plscmd','**','*'));dir(fullfile(root,'plsgui','**','*'))];
files=files(~[files.isdir]);
names=arrayfun(@(x)fullfile(x.folder,x.name),files,'UniformOutput',false);
names=sort(names);
entries=cell(numel(names),1);
for k=1:numel(names)
    relative=names{k}(numel(char(root))+2:end);
    entries{k}=[relative char(9) sha256_file(names{k})];
end
I=struct('root',char(root),'file_count',numel(names),'matlab_source_count',sum(endsWith(names,'.m')), ...
    'fingerprint',sha256_bytes(unicode2native(strjoin(entries,newline),'UTF-8')), ...
    'engine_sha256',sha256_file(fullfile(root,'plscmd','pls_analysis.m')));
end

function h = sha256_file(fp)
fid=fopen(fp,'rb'); assert(fid>=0,'Figure4:Unreadable','Cannot read identity input.');
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    bytes=fread(fid,1048576,'*uint8'); md.update(typecast(bytes,'int8'));
end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function h = sha256_bytes(bytes)
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(typecast(uint8(bytes(:)),'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end

function write_json(fp,value)
fid=fopen(fp,'w'); assert(fid>=0,'Figure4:Output','Cannot create requested report.');
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value));
end

function bsr = plsc_bootstrap_ratios(result, n_lv)
% Corrected frozen helper contract: compare_u IS the method-3 scaled-salience BSR.
if ~isfield(result,'boot_result') || ~isfield(result.boot_result,'compare_u') || ~isfield(result.boot_result,'u_se')
    error('PLSC:MissingBSR','Required saved compare_u/u_se fields missing; no fallback ratio.');
end
if ~isfield(result,'method') || result.method ~= 3, error('PLSC:Method','Method 3 required.'); end
u=double(result.u); s=double(result.s(:)); cu=double(result.boot_result.compare_u); se=double(result.boot_result.u_se);
if ~isreal(u) || ~isreal(cu) || ~isreal(se) || ~isequal(size(u),size(cu)) || ~isequal(size(u),size(se)) || numel(s)~=size(u,2)
    error('PLSC:BSRShape','Saved salience/singular-value/BSR/SE dimensions disagree.');
end
if n_lv<1 || n_lv~=fix(n_lv) || n_lv>size(u,2), error('PLSC:BSRLV','Invalid exported LV count.'); end
if any(~isfinite(u(:))) || any(~isfinite(cu(:))) || any(~isfinite(se(:))) || any(se(:)<=0) || any(~isfinite(s)) || any(s<=0)
    error('PLSC:BSRValues','Nonfinite salience/ratio or nonpositive SE/singular value.');
end
expected=bsxfun(@times,u,s')./se;
if isfield(result.boot_result,'zero_u_se') && ~isempty(result.boot_result.zero_u_se)
    z=double(result.boot_result.zero_u_se(:));
    if any(~isfinite(z)) || any(z~=fix(z)) || any(z<1) || any(z>numel(u)), error('PLSC:ZeroSE','Invalid zero-SE indices.'); end
    if any(se(z)~=1) || any(cu(z)~=0), error('PLSC:ZeroSE','Invalid toolbox zero-SE sentinel.'); end
    expected(z)=0;
end
if any(abs(cu(:)-expected(:))>1e-8+5e-6*abs(expected(:)))
    error('PLSC:BSRIdentity','compare_u differs from method-3 u*diag(s)/u_se.');
end
bsr=cu(:,1:n_lv);
end
