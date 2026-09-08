'use strict';
const D=JSON.parse(document.getElementById('report-data').textContent),$=id=>document.getElementById(id),NS='http://www.w3.org/2000/svg';
const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));

const tableRegistry=new Map();
const human=k=>String(k).replaceAll('_',' ');
const statusLabels={pass:'✓ Complete / supported',warn:'! Review',fail:'× Failed',not_evaluated:'— Not evaluated',info:'i Descriptive'};
function badge(state){const key=Object.hasOwn(statusLabels,state)?state:'info';return '<span class="status '+key+'">'+esc(statusLabels[key])+'</span>';}
function valueText(v){return v===null||v===undefined||v===''?'not recorded':Array.isArray(v)?(v.length?JSON.stringify(v):'none recorded'):typeof v==='object'?JSON.stringify(v):String(v);}

function displayCell(key,value){
 if(key==='state')return badge(value);
 const quality={high_confidence:'pass',moderate_confidence:'info',low_confidence:'warn',uncertain:'not_evaluated',rejected:'fail'};
 if(['quality_status','declared_quality_status'].includes(key)&&quality[value])return '<span class="status '+quality[value]+'">'+esc(human(value))+'</span>';
 if(key==='status'&&['complete','failed','not_evaluated'].includes(value))return '<span class="status '+({complete:'pass',failed:'fail',not_evaluated:'not_evaluated'}[value])+'">'+esc(human(value))+'</span>';
 return esc(valueText(value));
}

function downloadText(name,content,type='text/plain;charset=utf-8'){const url=URL.createObjectURL(new Blob([content],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),2000);}
function tsvText(rows,keys){const cell=v=>{let s=String(v??'');if(/^[=+@-]/.test(s)&&!/^[-+]?\d+(\.\d+)?$/.test(s))s="'"+s;return /[\t\r\n"]/.test(s)?'"'+s.replaceAll('"','""')+'"':s;};return [keys.map(cell).join('\t'),...rows.map(r=>keys.map(k=>cell(r[k])).join('\t'))].join('\n')+'\n';}
function table(id,rows,keys,onRow){
 const root=$(id);if(!root)return;
 const state={query:'',key:null,descending:false,page:0,size:25,rows,keys,onRow,printing:false};
 tableRegistry.set(id,state);
 root.innerHTML='<div class="table-tools"><label>Search this table<input type="search" aria-label="Search '+esc(human(id))+'"></label><button class="table-export">Export filtered TSV</button></div><div class="table-scroll"></div><div class="table-foot"><span aria-live="polite"></span><div><button class="prev">Previous</button><button class="next">Next</button></div></div>';
 const filtered=()=>{let data=rows.filter(r=>!state.query||keys.some(k=>String(r[k]??'').toLowerCase().includes(state.query)));
 if(state.key)data=data.slice().sort((a,b)=>{const av=a[state.key],bv=b[state.key],an=Number(av),bn=Number(bv),numeric=av!==''&&bv!==''&&av!=null&&bv!=null&&Number.isFinite(an)&&Number.isFinite(bn);return (numeric?an-bn:String(av??'').localeCompare(String(bv??''),undefined,{numeric:true}))*(state.descending?-1:1);});return data;};
 function draw(){const data=filtered(),pages=Math.max(1,Math.ceil(data.length/state.size));state.page=Math.min(state.page,pages-1);const shown=state.printing?data:data.slice(state.page*state.size,(state.page+1)*state.size);
 const container=root.querySelector('.table-scroll');
 container.innerHTML='<table><thead><tr>'+keys.map(k=>'<th scope="col" aria-sort="'+(state.key===k?(state.descending?'descending':'ascending'):'none')+'"><button data-key="'+esc(k)+'">'+esc(human(k))+(state.key===k?(state.descending?' ↓':' ↑'):' ↕')+'</button></th>').join('')+'</tr></thead><tbody>'+shown.map(r=>'<tr'+(onRow?' class="clickable" tabindex="0"':'')+'>'+keys.map(k=>'<td>'+displayCell(k,r[k])+'</td>').join('')+'</tr>').join('')+'</tbody></table>'+(data.length?'':'<p class="empty">No records for this selection. Missing analysis is not a negative biological result.</p>');
 container.querySelectorAll('thead button').forEach(b=>b.onclick=()=>{state.descending=state.key===b.dataset.key?!state.descending:false;state.key=b.dataset.key;state.page=0;draw();});
 if(onRow)container.querySelectorAll('tbody tr').forEach((el,i)=>{el.onclick=()=>onRow(shown[i]);el.onkeydown=e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();onRow(shown[i]);}};});
 root.querySelector('.table-foot span').textContent=data.length+' / '+rows.length+' rows · page '+(state.page+1)+' of '+pages;
 root.querySelector('.prev').disabled=state.page===0;root.querySelector('.next').disabled=state.page>=pages-1;}
 state.draw=draw;root.querySelector('input').oninput=e=>{state.query=e.target.value.toLowerCase();state.page=0;draw();};
 root.querySelector('.prev').onclick=()=>{state.page--;draw();};root.querySelector('.next').onclick=()=>{state.page++;draw();};
 root.querySelector('.table-export').onclick=()=>downloadText(id+'.tsv',tsvText(filtered(),keys),'text/tab-separated-values;charset=utf-8');draw();
}
function details(id,obj){$(id).innerHTML='<dl>'+Object.entries(obj).map(([k,v])=>'<dt>'+esc(human(k))+'</dt><dd>'+esc(valueText(v))+'</dd>').join('')+'</dl>';}
function options(id,values){$(id).innerHTML=values.map(x=>'<option value="'+esc(x)+'">'+esc(x)+'</option>').join('');}
function svg(tag,attrs,parent){const el=document.createElementNS(NS,tag);for(const [k,v] of Object.entries(attrs))el.setAttribute(k,v);parent.appendChild(el);return el;}
const unique=xs=>[...new Set(xs)].sort(),palette=['#087e83','#b34e32','#6255a1','#947125','#287245','#a74477','#485d88'];
const M=Object.fromEntries(D.metadata.map(r=>[r.isolate_id,r]));
$('run-label').textContent=D.provenance.dataset_kind+' · v'+D.provenance.framework_version+' · '+D.provenance.distance_engine+' · threshold '+D.provenance.threshold;
$('notice').textContent=D.provenance.dataset_kind.includes('synthetic')?'SYNTHETIC DEMONSTRATION — sequences are generated fixtures. Functional labels are illustrative, not detections or biological findings.':'Research evidence: review biological quality, annotation completion and epidemiological context. Shared units do not prove direct transmission.';
$('overview-cards').innerHTML=Object.entries(D.counts).filter(([k])=>!['cross_cluster_unit_links','typed_isolates'].includes(k)).map(([k,v])=>'<div class="card"><b>'+esc(k==='cross_cluster_pairs'&&!D.counts.typed_isolates?'—':v??'—')+'</b>'+esc(k.replaceAll('_',' '))+(k==='cross_cluster_pairs'&&!D.counts.typed_isolates?' · not evaluated':'')+'</div>').join('');
$('annotation-state').textContent='Annotation status: '+D.provenance.annotation_status+'. '+D.counts.typed_isolates+'/'+D.counts.isolates+' isolates have known chromosome clusters. Headline ARG counts describe eligible caller observations; stage completion is shown below.';
table('cross-table',D.crosslinks,['plasmid_unit','isolate_a','chrom_cluster_a','isolate_b','chrom_cluster_b']);
table('quality-table',D.plasmids,['plasmid_id','isolate_id','source_type','source_tool','circularity_status','declared_quality_status','quality_status','quality_warnings','duplicate_of'],r=>details('quality-detail',r));
table('biological-table',D.biological_quality||[],['plasmid_id','classification','replicon_type','mobility_class','graph_closure','read_breadth','mean_depth','chromosome_fraction','ambiguous_fraction','terminal_overlap_bp','quality_status','quality_warnings'],r=>details('quality-detail',r));
table('annotation-table',D.annotation_status||[],['plasmid_id','stage','status','engine','engine_version','database_version','cache','n_features']);
table('sensitivity-table',D.sensitivity,['threshold','linkage','n_units','n_singletons']);
$('safeguards').innerHTML=D.safeguards.map(s=>'<li>'+esc(s)+'</li>').join('');$('provenance').textContent=JSON.stringify(D.provenance,null,2);
const args=unique(D.features.filter(f=>f.headline_eligible==='true').map(f=>f.amr_gene));$('arg-filter').innerHTML+=args.map(a=>'<option>'+esc(a)+'</option>').join('');
let visibleEdges=[];
function network(){
 const query=$('search').value.toLowerCase(),arg=$('arg-filter').value,by=$('colour').value;
 const carriers=new Set(D.features.filter(f=>f.headline_eligible==='true'&&f.amr_gene===arg).map(f=>f.isolate_id));
 const matchingNodes=D.metadata.filter(r=>(!query||[r.isolate_id,r.location,r.organism].join(' ').toLowerCase().includes(query))&&(!arg||carriers.has(r.isolate_id)));
 const nodes=matchingNodes.slice(0,200);
 const present=new Set(nodes.map(r=>r.isolate_id)),groups=unique(D.metadata.map(r=>r[by]||'Unknown'));
 const colour=v=>palette[groups.indexOf(v||'Unknown')%palette.length];
 const cross=e=>M[e.source]?.chromosomal_cluster&&M[e.target]?.chromosomal_cluster&&M[e.source].chromosomal_cluster!==M[e.target].chromosomal_cluster;
 visibleEdges=D.edges.filter(e=>present.has(e.source)&&present.has(e.target)&&($('cross-only').value!=='cross'||cross(e))&&(!arg||D.edge_evidence.some(x=>x.source===e.source&&x.target===e.target&&x.shared_args.split(';').includes(arg))));
 const root=$('network');root.replaceChildren();const positions=Object.create(null);nodes.forEach((n,i)=>{const theta=2*Math.PI*i/Math.max(nodes.length,1)-Math.PI/2;positions[n.isolate_id]=[320+235*Math.cos(theta),210+155*Math.sin(theta)];});
 visibleEdges.forEach(e=>{const a=positions[e.source],b=positions[e.target],el=svg('line',{x1:a[0],y1:a[1],x2:b[0],y2:b[1],stroke:cross(e)?'#b34e32':'#9cadb1','stroke-width':2+Math.log2(Number(e.weight)),tabindex:0,role:'button','aria-label':e.source+' and '+e.target},root);svg('title',{},el).textContent=e.source+' ↔ '+e.target;el.onclick=()=>details('edge-detail',{...e,source_metadata:M[e.source],target_metadata:M[e.target],evidence:D.edge_evidence.filter(x=>x.source===e.source&&x.target===e.target),interpretation:'Candidate plasmid-sharing link; no direction or direct transmission established'});el.onkeydown=ev=>{if(ev.key==='Enter')el.onclick();};});
 nodes.forEach(n=>{const [x,y]=positions[n.isolate_id],el=svg('circle',{cx:x,cy:y,r:13,fill:colour(n[by]),stroke:'white','stroke-width':3,tabindex:0,role:'button','aria-label':n.isolate_id},root);const show=()=>details('edge-detail',{...n,plasmids:D.plasmids.filter(p=>p.isolate_id===n.isolate_id).map(p=>({id:p.plasmid_id,unit:p.plasmid_unit,quality:p.quality_status})),observation:'No indexed plasmid is not proof that the isolate is plasmid-free'});el.onclick=show;el.onkeydown=e=>{if(e.key==='Enter')show();};svg('text',{x,y:y+30,'text-anchor':'middle',fill:'#172f38'},root).textContent=n.isolate_id;});
 $('network-count').textContent=nodes.length+' visible isolates · '+visibleEdges.length+' visible sharing pairs'+(matchingNodes.length>200?' · display limited to the first 200 matches; narrow the search.':'');
 const typed=nodes.filter(n=>n.chromosomal_cluster).length,xc=visibleEdges.filter(cross).length;
 $('network-interpretation').textContent='Current network view: '+nodes.length+' isolates, '+typed+' with chromosome typing, '+visibleEdges.length+' sharing pairs and '+xc+' pairs across known clusters. '+(!typed?'Chromosome discordance is not evaluated in this view. ':'')+'These are candidate sharing links; direction and direct transmission remain unproven.';
 $('legend').innerHTML=groups.map(g=>'<span style="margin-right:16px;color:'+colour(g)+'">● '+esc(g)+'</span>').join('');
}
['search','arg-filter','colour','cross-only'].forEach(id=>{for(const event of ['input','change'])$(id).addEventListener(event,network);});
$('export-edges').onclick=()=>downloadText('visible_edges.tsv',tsvText(visibleEdges,['source','target','weight','shared_units']),'text/tab-separated-values;charset=utf-8');
const units=unique(D.plasmids.map(p=>p.plasmid_unit).filter(Boolean));options('unit',units);
const categories=unique(D.features.map(f=>f.functional_category||'unclassified'));$('category').innerHTML+=categories.map(c=>'<option>'+esc(c)+'</option>').join('');
const categoryColours={amr:'#b34e32',replication:'#087e83',mobility:'#6255a1',conjugation:'#6255a1',mobile_element:'#947125',metal_resistance:'#287245',metabolism:'#a74477',hypothetical:'#a6afb1'};
function featureDetail(f){const unit=f.plasmid_unit,neighbours=D.features.filter(x=>x.plasmid_id===f.plasmid_id&&x.source_sequence_id===f.source_sequence_id&&x.feature_id!==f.feature_id).sort((a,b)=>Math.abs(a.start-f.start)-Math.abs(b.start-f.start)).slice(0,4).map(x=>x.gene_symbol||x.feature_id);details('gene-detail',{...f,module_completeness:'unresolved',neighbouring_genes:neighbours,carrying_isolate:M[f.isolate_id],unit_isolates:unique(D.plasmids.filter(p=>unit&&p.plasmid_unit===unit).map(p=>p.isolate_id)),quality_warnings:D.plasmids.find(p=>p.plasmid_id===f.plasmid_id)?.quality_warnings});}
function cargo(){
 const pid=$('plasmid').value,unit=$('unit').value,category=$('category').value,query=$('gene-search').value.toLowerCase(),zoom=Number($('zoom').value);
 const rows=D.features.filter(f=>f.plasmid_id===pid&&(!category||(f.functional_category||'unclassified')===category)&&(!query||Object.values(f).join(' ').toLowerCase().includes(query)));
 const members=D.plasmids.filter(p=>unit&&p.plasmid_unit===unit),isos=unique(members.map(p=>p.isolate_id)),dates=isos.map(i=>M[i]?.date).filter(Boolean).sort();
 details('unit-detail',{plasmid_unit:unit||'no accepted units',representative_plasmid:members.map(p=>p.plasmid_id).sort()[0]||'',carrying_isolates:isos,organisms:unique(isos.map(i=>M[i]?.organism).filter(Boolean)),locations:unique(isos.map(i=>M[i]?.location).filter(Boolean)),date_range:dates.length?dates[0]+' to '+dates.at(-1):'not recorded',quality:'Inspect validation table; representative is lexicographically selected, not quality-ranked'});
 table('gene-table',rows,['gene_symbol','product_name','functional_category','start','end','strand','headline_eligible'],featureDetail);
 table('unit-table',D.unit_functions.filter(f=>f.plasmid_unit===unit),['feature_id','functional_category','prevalence_in_unit','n_plasmids','denominator_plasmids','n_evaluated_plasmids','annotation_coverage','n_isolates','completeness']);
 const root=$('track'),contigs=D.contigs[pid]||[],width=900*zoom,height=Math.max(100,contigs.length*100);root.replaceChildren();root.setAttribute('viewBox',`0 0 ${width} ${height}`);root.style.width=width+'px';root.style.height=height+'px';
 contigs.forEach((c,i)=>{const y=40+i*100;svg('text',{x:12,y:y-22,fill:'#52676e'},root).textContent=c.id+' · '+c.length+' bp';svg('line',{x1:20,y1:y,x2:width-25,y2:y,stroke:'#b8cacc','stroke-width':2},root);for(const f of rows.filter(f=>f.source_sequence_id===c.id)){const scale=(width-50)/c.length,x=20+(f.start-1)*scale,end=20+f.end*scale,tip=Math.min(9,(end-x)/2),left=f.strand==='-',points=f.strand==='.'?`${x},${y-9} ${end},${y-9} ${end},${y+9} ${x},${y+9}`:left?`${end},${y-9} ${x+tip},${y-9} ${x},${y} ${x+tip},${y+9} ${end},${y+9}`:`${x},${y-9} ${end-tip},${y-9} ${end},${y} ${end-tip},${y+9} ${x},${y+9}`;const arrow=svg('polygon',{points,fill:categoryColours[f.functional_category]||'#647789',tabindex:0,role:'button','aria-label':f.gene_symbol||f.feature_id},root);svg('title',{},arrow).textContent=(f.gene_symbol||f.feature_id)+' · '+f.product_name;arrow.onclick=()=>featureDetail(f);arrow.onkeydown=e=>{if(e.key==='Enter')featureDetail(f);};if(zoom>1||['amr','replication','mobility','conjugation'].includes(f.functional_category)){svg('text',{x,y:y+27,fill:'#172f38'},root).textContent=f.gene_symbol||f.feature_id;}}});
 if(!rows.length)svg('text',{x:20,y:height-12,fill:'#52676e'},root).textContent='No reported features for this selection; consult stage completion before interpretation.';
 const stages=(D.annotation_status||[]).filter(s=>s.plasmid_id===pid),done=stages.filter(s=>s.status==='complete').map(s=>s.stage),all=D.features.filter(f=>f.plasmid_id===pid),eligible=all.filter(f=>f.headline_eligible==='true');
 $('cargo-interpretation').textContent=pid?('Selected candidate '+pid+': '+all.length+' feature observations, '+eligible.length+' eligible ARG calls; '+rows.length+' features match the current filters. Completed stages: '+(done.join(', ')||'none recorded')+'. '+(!all.length?'No reported features do not establish biological absence. ':'')+'Product labels do not establish complete pathways or observed transfer.'):'No accepted candidate is available for functional inspection.';
 $('gene-detail').textContent='Select a feature to inspect provenance and metadata.';
}
function selectUnit(){options('plasmid',D.plasmids.filter(p=>p.plasmid_unit&&p.plasmid_unit===$('unit').value).map(p=>p.plasmid_id));cargo();}
$('unit').onchange=selectUnit;['plasmid','category','gene-search','zoom'].forEach(id=>{for(const event of ['input','change'])$(id).addEventListener(event,cargo);});

