
const A=D.dashboard,C=A.charts;
$('report-time').textContent='Generated '+new Date(D.report.generated_at).toLocaleString();
$('run-label').textContent=D.provenance.dataset_kind+' · analysis v'+D.provenance.framework_version+' / report v'+D.report.renderer_version+' · '+D.provenance.distance_engine+' · threshold '+D.provenance.threshold+' · '+D.provenance.status;
if(D.provenance.status!=='complete'){$('notice').classList.add('fail');$('notice').textContent='INCOMPLETE RUN — '+(D.provenance.error||D.provenance.status)+'. Available files are retained for diagnosis; missing modules are not negative results.';}
$('findings').innerHTML=A.findings.slice(0,6).map(f=>'<article class="finding '+esc(f.state)+'">'+badge(f.state)+'<h3>'+esc(f.title)+'</h3><p>'+esc(f.interpretation)+'</p><a href="#'+esc(f.section)+'">Review evidence →</a><small>Rule '+esc(f.rule_id)+' · v'+esc(f.rule_version)+'</small></article>').join('');
table('module-table',A.modules,['module','state','interpretation','evidence','next_step'],r=>{$(r.section).scrollIntoView({block:'start'});});
table('sample-table',A.samples,['isolate_id','state','accepted_candidates','rejected_candidates','total_length_bp','units','sharing_partners','amr_completed','annotation_denominator','eligible_arg_labels','chromosomal_cluster'],r=>selectIsolate(r.isolate_id));
table('metadata-table',D.metadata,unique(D.metadata.flatMap(r=>Object.keys(r))));
table('edge-table',D.edge_evidence,['source','target','plasmid_unit','minimum_distance','direct_threshold_support','shared_args']);
table('unit-membership',D.plasmids.filter(p=>p.plasmid_unit),['plasmid_unit','plasmid_id','isolate_id','length','n_contigs','quality_status','source_tool'],r=>selectCandidate(r.plasmid_id));
table('run-table',Object.entries(D.provenance).filter(([k,v])=>typeof v!=='object').map(([key,value])=>({key,value})),['key','value']);
const versions=new Map();for(const s of D.annotation_status||[]){if(s.engine)versions.set([s.engine,s.engine_version,s.database_sha256].join('|'),{engine:s.engine,version:s.engine_version,database_version:s.database_version,database_sha256:s.database_sha256});}
if(D.provenance.mash_version)versions.set('mash',{engine:'Mash',version:D.provenance.mash_version,database_version:'not applicable',database_sha256:'not applicable'});
table('versions-table',[...versions.values()],['engine','version','database_version','database_sha256']);
table('mobility-table',(D.biological_quality||[]).map(q=>({plasmid_id:q.plasmid_id,replicon_type:q.replicon_type,mobility_class:q.mobility_class,quality_status:q.quality_status})),['plasmid_id','replicon_type','mobility_class','quality_status']);
table('containment-table',D.containment||[],['small_plasmid_id','large_plasmid_id','small_isolate_id','large_isolate_id','small_length','large_length','length_ratio','pairwise_distance','max_distance']);
table('multilayer-table',D.multilayer_edges||[],['source','target','plasmid_unit','cluster_relation','minimum_distance','threshold_margin']);
table('typing-table',D.typing_crossreference||[],['plasmid_id','external_tool','mob_primary_cluster_id','ptu_assignment','predicted_host_range_overall_name','predicted_transmissibility_call','predicted_classification_call','plasmidfinder_inc_types','plsdb_nearest_accession','plsdb_nearest_distance','evidence_source']);
$('category-legend').innerHTML=Object.entries(categoryColours).map(([name,colour])=>'<span style="color:'+colour+'">● '+esc(human(name))+'</span>').join('');
$('status-palette').onchange=e=>{document.body.dataset.palette=e.target.value;};
$('print-report').onclick=()=>window.print();
window.addEventListener('beforeprint',()=>{document.querySelectorAll('main details').forEach(d=>{d.dataset.wasOpen=d.open?'yes':'no';d.open=true;});for(const t of tableRegistry.values()){t.printing=true;t.draw();}});
window.addEventListener('afterprint',()=>{document.querySelectorAll('main details').forEach(d=>{d.open=d.dataset.wasOpen==='yes';});for(const t of tableRegistry.values()){t.printing=false;t.draw();}});
if(!D.report.bundle){$('bundle-link').hidden=true;$('bundle-top').hidden=true;}
function exportSVG(root,name){const copy=root.cloneNode(true);copy.setAttribute('xmlns',NS);const originals=root.querySelectorAll('*'),clones=copy.querySelectorAll('*');originals.forEach((el,i)=>{const style=getComputedStyle(el);for(const key of ['fill','stroke','font-family','font-size'])clones[i].style.setProperty(key,style.getPropertyValue(key));});downloadText(name+'.svg',new XMLSerializer().serializeToString(copy),'image/svg+xml');}
function chartActions(container,root,rows,name){const actions=document.createElement('div');actions.className='chart-actions';const a=document.createElement('button');a.textContent='Download SVG';a.onclick=()=>exportSVG(root,name);const b=document.createElement('button');b.textContent='All chart data TSV';b.onclick=()=>downloadText(name+'.tsv',tsvText(rows,unique(rows.flatMap(r=>Object.keys(r)))),'text/tab-separated-values;charset=utf-8');actions.append(a,b);container.append(actions);}
function barChart(id,all,opts={}){
 const container=$(id),rows=all.slice(0,opts.limit||20),scaleMode=opts.scale||'linear',tf=v=>scaleMode==='log'?Math.log10(Number(v)+1):Number(v);
 container.replaceChildren();
 if(!rows.length){container.innerHTML='<p class="empty">No observations available. See module status for evaluation coverage.</p>';return;}
 const width=520,height=rows.length*33+35,linearMax=Math.max(1,...rows.map(r=>Number(r.denominator??r.value))),max=Math.max(1,...rows.map(r=>tf(r.denominator??r.value))),labelWidth=172,plotWidth=285;
 const root=svg('svg',{viewBox:'0 0 '+width+' '+height,role:'img','aria-label':opts.title||human(id),class:'chart'},container);
 rows.forEach((r,i)=>{const y=i*33+10,fill=opts.quality?({high_confidence:'var(--pass)',moderate_confidence:'var(--info)',low_confidence:'var(--warn)',uncertain:'var(--neutral)',rejected:'var(--fail)'}[r.label]||'var(--info)'):'#187c90';
 const text=svg('text',{x:0,y:y+14,fill:'#294859'},root);text.textContent=human(r.label).length>25?human(r.label).slice(0,23)+'…':human(r.label);svg('title',{},text).textContent=human(r.label);
 svg('rect',{x:labelWidth,y,width:plotWidth,height:20,rx:3,fill:'#edf2f5'},root);
 const mark=svg('rect',{x:labelWidth,y,width:plotWidth*tf(r.value)/max,height:20,rx:3,fill,tabindex:0,role:'button',class:'chart-mark','aria-label':r.label+': '+r.value+(r.denominator!==undefined?' of '+r.denominator:'')},root);
 const description=r.label+': '+r.value+(r.denominator!==undefined?' / '+r.denominator:'')+'. '+(opts.note||'Count from the full cohort.');
 svg('title',{},mark).textContent=description;mark.onclick=()=>{if(opts.onClick)opts.onClick(r);else $('chart-description-'+id).textContent=description;};mark.onkeydown=e=>{if(e.key==='Enter')mark.onclick();};
 svg('text',{x:labelWidth+plotWidth+7,y:y+14,fill:'#294859'},root).textContent=String(r.value)+(r.denominator!==undefined?'/'+r.denominator:'');});
 svg('text',{x:labelWidth,y:height-2,fill:'#526875'},root).textContent='0';svg('text',{x:labelWidth+plotWidth,y:height-2,'text-anchor':'end',fill:'#526875'},root).textContent=linearMax+' candidates / observations'+(scaleMode==='log'?' (log-scaled bar lengths)':'');
 const desc=document.createElement('p');desc.className='muted';desc.id='chart-description-'+id;desc.setAttribute('aria-live','polite');desc.textContent=all.length>rows.length?'Showing '+rows.length+' of '+all.length+' categories; export includes every category.':'Select a bar for its exact value.';container.append(desc);chartActions(container,root,all,id);
 if(opts.logToggle){const toggle=document.createElement('button');toggle.textContent=scaleMode==='log'?'Switch to linear scale':'Switch to log scale';toggle.title='Log view changes only the bar-length mapping for visual comparison across very different magnitudes; the printed value and TSV export always show the true linear count.';toggle.onclick=()=>barChart(id,all,{...opts,scale:scaleMode==='log'?'linear':'log'});container.querySelector('.chart-actions').append(toggle);}
}
barChart('chart-quality',C.quality,{quality:true,title:'Quality status of submitted candidates'});
barChart('chart-lengths',C.lengths,{title:'Lengths of accepted plasmid candidates'});
barChart('chart-annotation',C.annotation,{title:'Completed annotation stages by accepted candidate'});
barChart('chart-functions',C.functions,{title:'Observed functional category carriers',note:C.denominators.functions});
barChart('chart-metadata',C.metadata,{title:'Metadata coverage by isolate'});
barChart('chart-timeline',C.timeline,{title:'Isolate collection month counts',note:C.denominators.timeline});
barChart('chart-units',C.units,{title:'Candidate members per plasmid unit',logToggle:true,onClick:r=>{$('unit').value=r.label;selectUnit();$('cargo').scrollIntoView({block:'start'});}});
barChart('chart-drugs',C.drugs,{title:'Candidate carriers of eligible ARG drug classes',note:C.denominators.drugs});
function sensitivityChart(){
 const rows=D.sensitivity,container=$('chart-sensitivity');if(!rows.length){container.innerHTML='<p class="empty">No completed threshold sweep.</p>';return;}
 const width=540,height=245,left=50,right=510,top=20,bottom=205,xmax=Math.max(...rows.map(r=>Number(r.threshold)),.000001),ymax=Math.max(...rows.map(r=>Number(r.n_units)),1),root=svg('svg',{class:'chart',viewBox:'0 0 '+width+' '+height,role:'img','aria-label':'Unit counts by tested distance threshold'},container);
 const x=t=>left+Number(t)/xmax*(right-left),y=n=>bottom-Number(n)/ymax*(bottom-top);
 svg('path',{d:'M'+left+' '+top+'V'+bottom+'H'+right,fill:'none',stroke:'#9eb0bd'},root);
 for(const [key,color] of [['n_units','#137b87'],['n_singletons','#9b5e20']]){
 svg('polyline',{points:rows.map(r=>x(r.threshold)+','+y(r[key])).join(' '),fill:'none',stroke:color,'stroke-width':2},root);
 for(const r of rows){const mark=svg('circle',{cx:x(r.threshold),cy:y(r[key]),r:5,fill:color,tabindex:0,role:'button','aria-label':human(key)+' '+r[key]+' at threshold '+r.threshold},root);const caption=human(key)+': '+r[key]+' at distance '+r.threshold;svg('title',{},mark).textContent=caption;mark.onclick=()=>{$('sensitivity-caption').textContent=caption+'. This is sensitivity, not independent calibration.';};mark.onkeydown=e=>{if(e.key==='Enter')mark.onclick();};}}
 rows.forEach(r=>svg('text',{x:x(r.threshold),y:bottom+19,'text-anchor':'middle',fill:'#526875'},root).textContent=r.threshold);
 svg('text',{x:left-8,y:top+4,'text-anchor':'end',fill:'#526875'},root).textContent=ymax;
 svg('text',{x:left-8,y:bottom,'text-anchor':'end',fill:'#526875'},root).textContent='0';
 const p=document.createElement('p');p.id='sensitivity-caption';p.textContent='Teal: all units. Brown: singleton units. X: distance threshold; Y: unit count.';container.append(p);chartActions(container,root,rows,'threshold-sensitivity');
}
sensitivityChart();
function isolateTimeline(){
 const root=$('timeline');if(!root)return;root.replaceChildren();
 const width=960,height=260,left=40,right=920,base=200;
 root.setAttribute('viewBox','0 0 '+width+' '+height);
 const dated=D.metadata.filter(m=>m.date&&!isNaN(Date.parse(m.date))),undated=D.metadata.filter(m=>!(m.date&&!isNaN(Date.parse(m.date))));
 const datedIds=new Set(dated.map(m=>m.isolate_id));
 const shown=[],hidden=[];D.edge_evidence.forEach(e=>{(datedIds.has(e.source)&&datedIds.has(e.target)?shown:hidden).push(e);});
 if(!dated.length){svg('text',{x:20,y:100,fill:'#52676e'},root).textContent='No isolates with a usable collection date.';}
 else{
  const times=dated.map(m=>Date.parse(m.date)),xmin=Math.min(...times),xmax=Math.max(...times,xmin+1);
  const x=t=>left+(t-xmin)/(xmax-xmin)*(right-left),pos={};dated.forEach(m=>{pos[m.isolate_id]=x(Date.parse(m.date));});
  svg('line',{x1:left,y1:base,x2:right,y2:base,stroke:'#9eb0bd'},root);
  shown.forEach(e=>{const x1=pos[e.source],x2=pos[e.target],mid=(x1+x2)/2,arc=Math.min(90,Math.abs(x2-x1)/2+20),cross=crossCluster(e);
   const path=svg('path',{d:'M'+x1+' '+base+' Q '+mid+' '+(base-arc)+' '+x2+' '+base,fill:'none',stroke:cross?'#b34e32':'#9cadb1','stroke-width':2,tabindex:0,role:'button','aria-label':e.source+' to '+e.target},root);
   svg('title',{},path).textContent=e.source+' ↔ '+e.target+' · '+e.plasmid_unit+' · minimum distance '+e.minimum_distance;
   svg('text',{x:mid,y:base-arc-4,'text-anchor':'middle',fill:'#526875',style:'font-size:10px'},root).textContent=e.minimum_distance;});
  dated.forEach(m=>{const cx=pos[m.isolate_id],el=svg('circle',{cx,cy:base,r:6,fill:'#187c90',stroke:'white','stroke-width':2,tabindex:0,role:'button','aria-label':m.isolate_id},root);svg('title',{},el).textContent=m.isolate_id+' · '+m.date;el.onclick=()=>selectIsolate(m.isolate_id);el.onkeydown=ev=>{if(ev.key==='Enter')selectIsolate(m.isolate_id);};svg('text',{x:cx,y:base+18,'text-anchor':'middle',fill:'#172f38',style:'font-size:10px'},root).textContent=m.isolate_id;});
 }
 $('timeline-undated-isolates').textContent=undated.length?undated.length+' isolate(s) without a usable collection date: '+undated.map(m=>m.isolate_id).join(', ')+' (not omitted from the cohort; see Individual isolates).':'';
 table('timeline-undated-edges',hidden,['source','target','plasmid_unit','minimum_distance']);
}
isolateTimeline();
const V=D.distance_view;
$('distance-note').textContent='Distance engine: '+D.provenance.distance_engine+'; k='+D.provenance.k+'; sketch size '+(D.provenance.sketch_size??'not applicable')+'; linkage '+D.provenance.linkage+'. Distance is not measured sequence identity or transmission probability.';
$('heatmap-scope').textContent=V.ids.length?(V.limited?'Display limited to the first '+V.limit+' of '+V.total+' candidates. ':'All '+V.ids.length+' candidates shown. ')+'Dark teal indicates smaller distance; pale cells indicate larger distance, scaled to this displayed matrix. Use pair selectors for exact values; the complete matrix is in the output explorer.':'No completed distance matrix is available.';
options('pair-a',V.ids);options('pair-b',V.ids);if(V.ids.length>1)$('pair-b').selectedIndex=1;
const canvas=$('distance-map'),ctx=canvas.getContext('2d'),n=V.ids.length;
if(n){const max=Math.max(...V.matrix.flat().map(Number),.000001),cell=canvas.width/n;
 for(let i=0;i<n;i++)for(let j=0;j<n;j++){const t=Number(V.matrix[i][j])/max;ctx.fillStyle='rgb('+Math.round(12+225*t)+','+Math.round(93+152*t)+','+Math.round(112+137*t)+')';ctx.fillRect(j*cell,i*cell,Math.ceil(cell),Math.ceil(cell));}
 canvas.onclick=e=>{const box=canvas.getBoundingClientRect(),j=Math.min(n-1,Math.floor((e.clientX-box.left)/box.width*n)),i=Math.min(n-1,Math.floor((e.clientY-box.top)/box.height*n));$('pair-a').value=V.ids[i];$('pair-b').value=V.ids[j];showPair();};
}else canvas.hidden=true;
function showPair(){const i=V.ids.indexOf($('pair-a').value),j=V.ids.indexOf($('pair-b').value);$('pair-detail').textContent=i<0||j<0?'Pair distance not evaluated.':V.ids[i]+' ↔ '+V.ids[j]+': distance '+V.matrix[i][j]+'. '+(i===j?'This is a self-comparison.':Number(V.matrix[i][j])<=Number(D.provenance.threshold)?'Meets the current pairwise distance cutoff.':'Exceeds the current pairwise distance cutoff.')+' Pairwise support does not establish direct transmission; unit membership also depends on linkage and the cohort.';}
$('pair-a').onchange=showPair;$('pair-b').onchange=showPair;showPair();
function selectIsolate(iso){$('isolate-select').value=iso;renderIsolate();$('individual').scrollIntoView({block:'start'});}
function selectCandidate(pid){const p=D.plasmids.find(p=>p.plasmid_id===pid);if(!p?.plasmid_unit)return;$('unit').value=p.plasmid_unit;selectUnit();$('plasmid').value=pid;cargo();$('cargo').scrollIntoView({block:'start'});}
options('isolate-select',A.samples.map(r=>r.isolate_id));
function renderIsolate(){
 const iso=$('isolate-select').value,s=A.samples.find(r=>r.isolate_id===iso);
 if(!s){$('isolate-interpretation').textContent='No isolate metadata is available.';return;}
 $('isolate-interpretation').textContent=s.interpretation+' This summary concerns the selected isolate; cohort charts remain unchanged.';
 details('isolate-detail',{...M[iso],submitted_candidates:s.submitted_candidates,accepted_candidates:s.accepted_candidates,rejected_candidates:s.rejected_candidates,total_length_bp:s.total_length_bp,eligible_arg_labels:s.eligible_arg_labels||'none observed; check caller completion'});
 const ps=D.plasmids.filter(p=>p.isolate_id===iso),ids=new Set(ps.map(p=>p.plasmid_id));
 table('isolate-plasmids',ps,['plasmid_id','plasmid_unit','length','n_contigs','quality_status','circularity_status','source_tool','quality_warnings'],p=>selectCandidate(p.plasmid_id));
 table('isolate-stages',(D.annotation_status||[]).filter(r=>ids.has(r.plasmid_id)),['plasmid_id','stage','status','engine','engine_version','database_version','cache','n_features']);
 table('isolate-edges',D.edge_evidence.filter(e=>e.source===iso||e.target===iso),['source','target','plasmid_unit','minimum_distance','direct_threshold_support','shared_args']);
 const files=D.artifacts.filter(f=>ps.some(p=>f.path.startsWith('annotation/'+p.plasmid_id+'/')||f.path==='plasmids/'+p.plasmid_id+'.fasta'));
 table('isolate-files',files,['path','module','description','size_bytes'],f=>{inspectFile(f);$('downloads').scrollIntoView({block:'start'});});
}
$('isolate-select').onchange=renderIsolate;
$('export-sample').onclick=()=>{const iso=$('isolate-select').value;downloadText((iso||'isolate')+'.summary.json',JSON.stringify({summary:A.samples.find(r=>r.isolate_id===iso),metadata:M[iso],candidates:D.plasmids.filter(p=>p.isolate_id===iso)},null,2),'application/json');};
if(D.calibration){const c=D.calibration.record;
 $('calibration-interpretation').textContent='Status: '+c.status+'. Target: '+c.target+'. Training selected '+(c.selected_threshold??'no threshold')+'; current run uses '+D.provenance.threshold+'. '+c.limitation;
 details('calibration-detail',c);table('calibration-table',D.calibration.metrics,['split','threshold','tp','fp','fn','tn','precision','recall','specificity','f1','balanced_accuracy']);
}else{$('calibration-interpretation').textContent='Not evaluated in this report: no matching calibration attachment. Sensitivity is not a replacement for independent labels.';$('calibration-detail').textContent='Attach an existing matching result with plasmicord report --results RESULTS --calibration CALIBRATION. This does not recompute distances or change the threshold.';}
const bytes=n=>n==null?'final report file':n<1024?n+' B':n<1048576?(n/1024).toFixed(1)+' KiB':n<1073741824?(n/1048576).toFixed(1)+' MiB':(n/1073741824).toFixed(2)+' GiB';
let activeFile=null,expandFiles=false;
$('file-module').innerHTML+=unique(D.artifacts.map(f=>f.module)).map(m=>'<option>'+esc(m)+'</option>').join('');
function inspectFile(file){
 activeFile=file.path;$('file-title').textContent=file.path;
 $('file-metadata').textContent=file.description+'\nModule: '+file.module+' · Size: '+bytes(file.size_bytes)+' · Modified: '+(file.modified_utc||'not recorded')+'\nSHA-256: '+(file.sha256||'See output_manifest.json for finalized report hashes. The manifest and ZIP exclude their own hashes.');
 $('file-metadata').style.whiteSpace='pre-wrap';
 const actions=$('file-actions');actions.replaceChildren();
 const download=document.createElement('a');download.href=file.href;download.download=file.name;download.className='button';download.textContent='Download original';
 const open=document.createElement('a');open.href=file.href;open.target='_blank';open.rel='noopener noreferrer';open.className='button';open.textContent='Open original';
 actions.append(download,open);
 $('file-preview').textContent=file.preview!==null&&file.preview!==undefined?file.preview+(file.preview_truncated?'\n\n[Preview truncated. Download the complete original file.]':''):'No embedded text preview. Binary, large, or final report files remain available through the links above.';
 document.querySelectorAll('#file-tree button').forEach(b=>b.classList.toggle('active',b.dataset.path===activeFile));
}
function renderTree(){
 const query=$('file-search').value.toLowerCase(),module=$('file-module').value,files=D.artifacts.filter(f=>(!module||f.module===module)&&(!query||[f.path,f.description,f.module].join(' ').toLowerCase().includes(query)));
 const tree={directories:new Map(),files:[]};
 for(const file of files){const parts=file.path.split('/');let node=tree;for(const folder of parts.slice(0,-1)){if(!node.directories.has(folder))node.directories.set(folder,{directories:new Map(),files:[]});node=node.directories.get(folder);}node.files.push(file);}
 function render(node){const list=document.createElement('ul');for(const [name,child]of [...node.directories.entries()].sort((a,b)=>a[0].localeCompare(b[0]))){const li=document.createElement('li'),details=document.createElement('details'),label=document.createElement('summary');label.textContent=name+'/';details.open=expandFiles||!!query||!!module;details.append(label,render(child));li.append(details);list.append(li);}
 for(const f of node.files){const li=document.createElement('li'),button=document.createElement('button');button.textContent=f.name;button.dataset.path=f.path;button.classList.toggle('active',f.path===activeFile);button.onclick=()=>inspectFile(f);const small=document.createElement('span');small.className='file-size';small.textContent=bytes(f.size_bytes);li.append(button,small);list.append(li);}return list;}
 $('file-tree').replaceChildren(render(tree));$('file-count').textContent=files.length+' of '+D.artifacts.length+' files shown · '+bytes(files.reduce((sum,f)=>sum+(f.size_bytes||0),0))+' bytes measured before packaging. Final metadata and checksums are in output_manifest.json. All paths are relative to this report.';
}
$('file-search').oninput=renderTree;$('file-module').onchange=renderTree;
$('expand-tree').onclick=()=>{expandFiles=true;document.querySelectorAll('#file-tree details').forEach(d=>d.open=true);};
$('collapse-tree').onclick=()=>{expandFiles=false;document.querySelectorAll('#file-tree details').forEach(d=>d.open=false);};
const observer=new IntersectionObserver(entries=>{for(const entry of entries)if(entry.isIntersecting){document.querySelectorAll('.sidebar nav a').forEach(a=>a.classList.toggle('active',a.getAttribute('href')==='#'+entry.target.id));}},{rootMargin:'-5% 0px -70% 0px'});
document.querySelectorAll('main section').forEach(s=>observer.observe(s));
network();selectUnit();renderTree();renderIsolate();
