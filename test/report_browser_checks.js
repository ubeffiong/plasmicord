(() => {
  const assert = (ok, message) => { if (!ok) throw new Error(message); };
  const data = JSON.parse(document.getElementById('report-data').textContent);
  const node = id => document.getElementById(id);
  const change = (id, value, type='change') => { node(id).value=value; node(id).dispatchEvent(new Event(type,{bubbles:true})); };
  assert(node('file-tree').querySelectorAll('button').length === data.artifacts.length, 'Every inventoried output must be in the tree');
  assert(node('sample-table').querySelectorAll('tbody tr').length === Math.min(25,data.metadata.length), 'Sample pagination count');
  if(data.metadata.length > 25) {
    node('sample-table').querySelector('.next').click();
    assert(node('sample-table').querySelectorAll('tbody tr').length === data.metadata.length-25, 'Second page row count');
    node('sample-table').querySelector('.prev').click();
  }
  const search = node('sample-table').querySelector('input');
  if(data.metadata.length) {
    search.value=data.metadata[0].isolate_id; search.dispatchEvent(new Event('input',{bubbles:true}));
    assert([...node('sample-table').querySelectorAll('tbody tr')].every(r=>r.textContent.includes(data.metadata[0].isolate_id)), 'Table filtering');
    node('sample-table').querySelector('tbody tr').click();
    assert(node('isolate-select').value===data.metadata[0].isolate_id, 'Row opens correct isolate after filtering');
    search.value='';search.dispatchEvent(new Event('input',{bubbles:true}));
  }
  change('search','__no_such_isolate__','input');
  assert(node('network-count').textContent.startsWith('0 visible isolates'), 'Network empty state');
  assert(node('network-interpretation').textContent.includes('not evaluated'), 'Unknown view must not imply concordance');
  change('search','','input');
  change('status-palette','accessible');
  assert(document.body.dataset.palette==='accessible','Accessible palette');
  change('status-palette','standard');
  if(data.distance_view.ids.length>1) {
    change('pair-a',data.distance_view.ids[0]);change('pair-b',data.distance_view.ids[1]);
    assert(node('pair-detail').textContent.includes(String(data.distance_view.matrix[0][1])), 'Exact pairwise distance');
  }
  const nested=data.artifacts.find(f=>f.path.includes('/')&&f.preview!==null) || data.artifacts.find(f=>f.preview!==null);
  if(nested) {
    change('file-search',nested.path,'input');
    [...node('file-tree').querySelectorAll('button')].find(b=>b.dataset.path===nested.path).click();
    assert(node('file-title').textContent===nested.path,'Nested output selection');
    assert(node('file-preview').textContent.startsWith(nested.preview),'Text preview matches stored evidence');
    assert(node('file-actions').querySelector('a[download]').getAttribute('href')===nested.href,'Download URL');
  }
  change('file-search','','input');
  if(data.features.length) {
    const feature=data.features.find(f=>f.headline_eligible==='true')||data.features[0];
    change('unit',feature.plasmid_unit);change('plasmid',feature.plasmid_id);
    change('category',feature.functional_category);
    assert(node('cargo-interpretation').textContent.includes(feature.plasmid_id),'Selected candidate interpretation');
    const query=node('gene-table').querySelector('input');query.value=feature.gene_symbol||feature.product_name;query.dispatchEvent(new Event('input',{bubbles:true}));
    node('gene-table').querySelector('tbody tr').click();
    assert(node('gene-detail').textContent.includes(feature.annotation_engine),'Gene details retain caller provenance');
    if(node('track-view')) {
      assert(!node('track').hidden && node('circular-track').hidden,'Linear track view is the default');
      change('track-view','circular');
      assert(node('track').hidden && !node('circular-track').hidden,'Circular view toggle swaps visible tracks');
      assert(node('circular-track').querySelectorAll('circle').length>0,'Circular view draws a backbone ring per contig');
      const wedge=node('circular-track').querySelector('path');
      if(wedge) {
        const wedgeLabel=wedge.getAttribute('aria-label');
        wedge.dispatchEvent(new Event('click',{bubbles:true}));
        assert(node('gene-detail').textContent.includes(wedgeLabel),'Clicking a circular-map wedge opens that exact feature\'s detail, not a placeholder or a different feature');
      }
      change('track-view','linear');
      assert(!node('track').hidden && node('circular-track').hidden,'Toggling back restores the linear view');
    }
  }
  if(data.calibration) assert(node('calibration-table').textContent.includes('validation'), 'Holdout metrics surfaced');
  const charts=document.querySelectorAll('.chart').length;
  assert(charts>=4,'Chart modules rendered');
  if(node('size-by')) {
    const before=node('legend').innerHTML;
    assert(!before.includes('Node size'),'Node-size legend absent by default');
    change('size-by','degree');
    assert(node('legend').innerHTML.includes('Node size'),'Node-size legend appears once a size mode is active');
    change('size-by','none');
    assert(!node('legend').innerHTML.includes('Node size'),'Node-size legend removed when reset to none');
  }
  const firstEdgeLine=node('network').querySelector('line');
  if(firstEdgeLine) {
    firstEdgeLine.dispatchEvent(new Event('click',{bubbles:true}));
    assert(node('edge-unit').options.length>=1,'Selecting an edge populates its shared-unit selector');
    assert(node('edge-tracks-svg').children.length>0,'Selecting an edge renders its gene-track comparison');
    assert(node('edge-tracks-caption').textContent.includes('not sequence alignment'),'Gene-track comparison carries its non-alignment caveat');
    if(node('edge-unit').options.length>1) {
      const before=node('edge-tracks-svg').innerHTML;
      change('edge-unit',node('edge-unit').options[1].value);
      assert(node('edge-tracks-svg').innerHTML!==before,'Switching the shared-unit selector redraws the comparison');
    }
  }
  const datedCount=data.metadata.filter(m=>m.date&&!isNaN(Date.parse(m.date))).length;
  assert(node('timeline').querySelectorAll('circle').length===datedCount,'Timeline plots exactly the dated isolates');
  const hasUndated=data.metadata.length>datedCount, undatedTextShown=node('timeline-undated-isolates').textContent.length>0;
  assert(undatedTextShown===hasUndated,'Undated isolates are listed, not silently dropped');
  if(node('chart-units').querySelector('.chart-actions')) {
    const findToggle=()=>[...node('chart-units').querySelectorAll('.chart-actions button')].find(b=>/log|linear/i.test(b.textContent));
    const toggle=findToggle();
    assert(!!toggle,'Log/linear scale toggle present on the unit-size chart');
    const before=toggle.textContent;
    toggle.click(); // barChart() fully redraws and replaces the button; re-query, don't reuse the old reference
    const after=findToggle();
    assert(!!after&&after.textContent!==before,'Toggle switches its own label after a click');
  }
  window.dispatchEvent(new Event('beforeprint'));
  assert(node('sample-table').querySelectorAll('tbody tr').length===data.metadata.length,'Print includes all sample rows');
  window.dispatchEvent(new Event('afterprint'));
  assert(node('sample-table').querySelectorAll('tbody tr').length===Math.min(25,data.metadata.length),'Pagination restored after print');
  return {passed:true,isolates:data.metadata.length,features:data.features.length,files:data.artifacts.length,charts,calibration:!!data.calibration};
})()
