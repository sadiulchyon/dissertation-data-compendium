// Ch4 January low-flow site extraction (browser-injected helpers in the GEE Code Editor tab).
// Saved 2026-10-22. Paste into the DevTools console of a logged-in Code Editor tab, then call
// window.__runSites([[sid,W,S,E,N],...]) (list from scripts/windows.py). Results: window.__res[sid]; bins download as SITE_scenes_v2.bin.

window.__runSites = async function(sites){
  for(const [sid,W,S,E,N] of sites){
    window.__res[sid]='running';
    try{ window.__res[sid]=await window.__site(sid,W,S,E,N); }catch(e){ window.__res[sid]='ERROR '+e; }
    await new Promise(r=>setTimeout(r,1500));
  }
};

window.__site = async function(sid, W,S,E,N){
  const reg = ee.Geometry.Rectangle([W,S,E,N], null, false);
  const COL={L5:'LANDSAT/LT05/C02/T1_L2',L7:'LANDSAT/LE07/C02/T1_L2',L8:'LANDSAT/LC08/C02/T1_L2',L9:'LANDSAT/LC09/C02/T1_L2'};
  const BLUE={L5:'SR_B1',L7:'SR_B1',L8:'SR_B2',L9:'SR_B2'};
  const BN={L5:['SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'],L7:['SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'],L8:['SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'],L9:['SR_B3','SR_B4','SR_B5','SR_B6','SR_B7']};
  const clearOf=im=>{const qa=im.select('QA_PIXEL'); return qa.bitwiseAnd(1<<0).eq(0).and(qa.bitwiseAnd(1<<1).eq(0)).and(qa.bitwiseAnd(1<<3).eq(0)).and(qa.bitwiseAnd(1<<4).eq(0)).and(qa.bitwiseAnd(1<<2).eq(0));};
  function prep(sen){
    return ee.ImageCollection(COL[sen]).filterBounds(reg).filter(ee.Filter.calendarRange(11,4,'month')).map(function(im){
      const c=clearOf(im);
      const f=c.unmask(0).rename('c').reduceRegion({reducer:ee.Reducer.mean(),geometry:reg,scale:60,maxPixels:1e8}).get('c');
      const b=im.select(BLUE[sen]).multiply(0.0000275).add(-0.2).updateMask(c).rename('b').reduceRegion({reducer:ee.Reducer.median(),geometry:reg,scale:60,maxPixels:1e8}).get('b');
      return ee.Feature(null,{id:im.get('system:index'),sensor:sen,date:im.date().format('YYYY-MM-dd'),clear:f,blue:b});
    });
  }
  const all=prep('L5').merge(prep('L7')).merge(prep('L8')).merge(prep('L9'));
  const rows=await new Promise((res,rej)=>all.filter(ee.Filter.gte('clear',0.6)).reduceColumns(ee.Reducer.toList(5),['sensor','date','clear','blue','id']).get('list').evaluate((v,e)=>e?rej(e):res(v)));
  const Sx={}; for(const [sen,date,clear,blue,id] of rows){ const y=+date.slice(0,4),m=+date.slice(5,7); const s=m>=11?y+1:y; if(blue===null) continue; (Sx[s]=Sx[s]||[]).push({sen,date,clear:+clear.toFixed(3),blue:+blue.toFixed(4),id}); }
  const picks={}, log=[];
  for(let s=1988;s<=2025;s++){
    const L=Sx[s]||[]; const ref=new Date(s+'-01-20'); L.forEach(o=>o.dd=Math.round((new Date(o.date)-ref)/864e5));
    const core=L.filter(o=>o.dd>=-36&&o.dd<=39);
    const best=(arr,th)=>{const c=arr.filter(o=>o.clear>=th); c.sort((a,b)=>a.blue-b.blue); return c[0];};
    // prefer non-L7 when equally good: L7 SLC-off after 2003 has gaps (clear already reflects that)
    const p=best(core,0.9)||best(core,0.8)||best(L,0.9)||best(L,0.8)||best(core,0.6)||best(L,0.6);
    if(p){picks[s]=p; log.push(`${s} ${p.sen} ${p.date} clear=${p.clear} blue=${p.blue} (n core ${core.length})`);} else log.push(`${s} NONE`);
  }
  const gsw=ee.Image('JRC/GSW1_4/GlobalSurfaceWater');
  async function grab(img){ const url=await new Promise((res,rej)=>img.getDownloadURL({region:reg,crs:'EPSG:32646',scale:30,format:'NPY'},(u,e)=>e?rej(e):res(u))); const r=await fetch(url); if(r.status!==200) throw new Error('status '+r.status+' '+(await r.text()).slice(0,200)); return new Uint8Array(await r.arrayBuffer()); }
  const seasons=Object.keys(picks).map(Number).sort();
  const bufs=[await grab(gsw.select(['occurrence','max_extent','seasonality']).unmask(0).toInt16())];
  for(let i=0;i<seasons.length;i+=6){
    const part=await Promise.all(seasons.slice(i,i+6).map(async s=>{
      const p=picks[s]; const im=ee.Image(COL[p.sen]+'/'+p.id);
      const sr=im.select(BN[p.sen],['GREEN','RED','NIR','SWIR1','SWIR2']).multiply(0.0000275).add(-0.2).multiply(10000);
      return grab(sr.addBands(clearOf(im).rename('CLEAR')).unmask(-32768).toInt16());
    }));
    bufs.push(...part);
  }
  const meta=seasons.map(s=>({season:s,...picks[s]}));
  const manifest=JSON.stringify({site:sid,region:[W,S,E,N],crs:'EPSG:32646',scale:30,selection:'clear>=0.9 in 15 Dec-28 Feb, lowest median blue SR (haze); fallbacks clear>=0.8, Nov-Apr',bands:['GREEN','RED','NIR','SWIR1','SWIR2','CLEAR'],gsw_bands:['occurrence','max_extent','seasonality'],scenes:meta,lens:bufs.map(b=>b.length)});
  const head=new TextEncoder().encode(manifest+'\n'); const total=head.length+bufs.reduce((a,b)=>a+b.length,0);
  const out=new Uint8Array(total); out.set(head,0); let o=head.length; for(const b of bufs){out.set(b,o); o+=b.length;}
  window.__blob=new Blob([out]); window.__blobName=sid+'_scenes_v2.bin'; window.__save();
  return log.join('\n')+'\nsaved '+window.__blobName+' '+total;
};

window.__pull = async function(site, W,S,E,N){
  const reg = ee.Geometry.Rectangle([W,S,E,N], null, false);
  const bands=['GREEN','RED','NIR','SWIR1','SWIR2','NDVI','MNDWI','NOBS'];
  const occ = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence').unmask(0).rename('OCC');
  const years=[]; for(let y=1988;y<=2025;y++) years.push(y);
  const bufs=[];
  for (let i=0;i<years.length;i+=8){
    const part = await Promise.all(years.slice(i,i+8).map(async y=>{
      const img = ee.Image(ROOT+'ls_dry/ls_dry_'+y).select(bands).addBands(occ).toInt16();
      const url = await new Promise((res,rej)=>img.getDownloadURL({region: reg, scale: 30, crs:'EPSG:32646', format:'NPY'}, (u,e)=>e?rej(e):res(u)));
      const r = await fetch(url); if(r.status!==200) throw new Error('status '+r.status+' year '+y);
      return new Uint8Array(await r.arrayBuffer());
    }));
    bufs.push(...part);
  }
  const manifest = JSON.stringify({site, region:[W,S,E,N], scale:30, crs:'EPSG:32646', years, bands:[...bands,'OCC'], note:'reflectance and indices x1e4 as stored; masked = -32768', lens: bufs.map(b=>b.length)});
  const head = new TextEncoder().encode(manifest+'\n');
  const total = head.length + bufs.reduce((s,b)=>s+b.length,0);
  const all = new Uint8Array(total); all.set(head,0); let o=head.length; for(const b of bufs){ all.set(b,o); o+=b.length; }
  window.__blob = new Blob([all]); window.__blobName = site+'_cube.bin';
  return total;
};

window.__save = function(){ const a=document.createElement('a'); a.href=URL.createObjectURL(window.__blob); a.download=window.__blobName; document.body.appendChild(a); a.click(); a.remove(); return window.__blobName; };

window.__merit = async function(list){
  const bufs=[], sites=[];
  for(const [sid,W,S,E,N] of list){
    const reg=ee.Geometry.Rectangle([W,S,E,N],null,false);
    const img=ee.Image('MERIT/Hydro/v1_0_1').select(['upa','viswth']).unmask(0).toFloat();
    const url=await new Promise((res,rej)=>img.getDownloadURL({region:reg,crs:'EPSG:32646',scale:30,format:'NPY'},(u,e)=>e?rej(e):res(u)));
    const r=await fetch(url); if(r.status!==200) throw new Error(await r.text());
    bufs.push(new Uint8Array(await r.arrayBuffer())); sites.push(sid);
  }
  const head=new TextEncoder().encode(JSON.stringify({sites,lens:bufs.map(b=>b.length),bands:['upa','viswth'],note:'MERIT Hydro v1.0.1 resampled to the site 30 m grid (nearest)'})+'\n');
  const tot=head.length+bufs.reduce((a,b)=>a+b.length,0); const out=new Uint8Array(tot); out.set(head,0); let o=head.length; for(const b of bufs){out.set(b,o);o+=b.length;}
  const a=document.createElement('a'); a.href=URL.createObjectURL(new Blob([out])); a.download='merit_'+sites.join('_')+'.bin'; document.body.appendChild(a); a.click(); a.remove();
  return tot;
};

// Added 25 Oct 2026. Same-season Landsat 8 / Landsat 7 pairs (1-16 days apart, seasons 2014-2021) for the sensor check:
// window.__runPairs([[sid,W,S,E,N],...]) -> SITE_pair.bin (read by scripts/pair_test.py and scripts/pair_test_harm.py).
window.__pair = async function(sid,W,S,E,N){
  const reg = ee.Geometry.Rectangle([W,S,E,N], null, false);
  const COL={L7:'LANDSAT/LE07/C02/T1_L2',L8:'LANDSAT/LC08/C02/T1_L2'};
  const BN={L7:['SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'],L8:['SR_B3','SR_B4','SR_B5','SR_B6','SR_B7']};
  const clearOf=im=>{const qa=im.select('QA_PIXEL'); return qa.bitwiseAnd(1<<0).eq(0).and(qa.bitwiseAnd(1<<1).eq(0)).and(qa.bitwiseAnd(1<<3).eq(0)).and(qa.bitwiseAnd(1<<4).eq(0)).and(qa.bitwiseAnd(1<<2).eq(0));};
  function prep(sen){ return ee.ImageCollection(COL[sen]).filterBounds(reg).filterDate('2012-11-01','2021-04-30').filter(ee.Filter.calendarRange(11,3,'month')).map(function(im){
      const f=clearOf(im).unmask(0).rename('c').reduceRegion({reducer:ee.Reducer.mean(),geometry:reg,scale:60,maxPixels:1e8}).get('c');
      return ee.Feature(null,{id:im.get('system:index'),sensor:sen,date:im.date().format('YYYY-MM-dd'),clear:f}); }); }
  const rows=await new Promise((res,rej)=>prep('L7').merge(prep('L8')).filter(ee.Filter.gte('clear',0.5)).reduceColumns(ee.Reducer.toList(4),['sensor','date','clear','id']).get('list').evaluate((v,e)=>e?rej(e):res(v)));
  const L7=rows.filter(r=>r[0]=='L7'), L8=rows.filter(r=>r[0]=='L8'&&r[2]>=0.85);
  const season=d=>{const y=+d.slice(0,4),m=+d.slice(5,7);return m>=11?y+1:y;};
  const best={};
  for(const a of L8) for(const b of L7){ const dd=Math.abs((new Date(a[1])-new Date(b[1]))/864e5); if(dd>16||season(a[1])!=season(b[1])) continue;
    const s=season(a[1]); const sc=dd/16 - b[2]; if(!best[s]||sc<best[s].sc) best[s]={sc,a,b,dd}; }
  const seasons=Object.keys(best).map(Number).sort();
  async function grab(img){ const url=await new Promise((res,rej)=>img.getDownloadURL({region:reg,crs:'EPSG:32646',scale:30,format:'NPY'},(u,e)=>e?rej(e):res(u))); const r=await fetch(url); if(r.status!==200) throw new Error('status '+r.status); return new Uint8Array(await r.arrayBuffer()); }
  const meta=[], bufs=[];
  for(const s of seasons){ for(const r of [best[s].a,best[s].b]){ const sen=r[0]; const im=ee.Image(COL[sen]+'/'+r[3]);
      const sr=im.select(BN[sen],['GREEN','RED','NIR','SWIR1','SWIR2']).multiply(0.0000275).add(-0.2).multiply(10000);
      bufs.push(await grab(sr.addBands(clearOf(im).rename('CLEAR')).unmask(-32768).toInt16())); meta.push({season:s,sen,date:r[1],clear:r[2],id:r[3]}); } }
  const head=new TextEncoder().encode(JSON.stringify({site:sid,region:[W,S,E,N],bands:['GREEN','RED','NIR','SWIR1','SWIR2','CLEAR'],scenes:meta,lens:bufs.map(b=>b.length)})+'\n');
  const tot=head.length+bufs.reduce((a,b)=>a+b.length,0); const out=new Uint8Array(tot); out.set(head,0); let o=head.length; for(const b of bufs){out.set(b,o);o+=b.length;}
  window.__blob=new Blob([out]); window.__blobName=sid+'_pair.bin'; window.__save();
  return seasons.map(s=>`${s} L8 ${best[s].a[1]} L7 ${best[s].b[1]} (${best[s].b[2].toFixed(2)}) dd${best[s].dd}`).join('\n');
};
window.__runPairs = async function(list){ for(const [sid,W,S,E,N] of list){ window.__res[sid+'p']='running'; try{window.__res[sid+'p']=await window.__pair(sid,W,S,E,N);}catch(e){window.__res[sid+'p']='ERROR '+e;} await new Promise(r=>setTimeout(r,1500)); } };
// Sentinel-2 screening strip: Jan-Mar median true colour for the given years, saved as s2_NAME.png (vertical filmstrip).
window.__s2strip = async function(name,lon,lat,km,yrs){
  const reg=ee.Geometry.Point([lon,lat]).buffer(km*500).bounds();
  const ims=yrs.map(y=>ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(reg).filterDate(y+'-01-01',y+'-03-31').filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',15)).median().visualize({bands:['B4','B3','B2'],min:300,max:2200}));
  const url=await new Promise((res,rej)=>ee.ImageCollection(ims).getFilmstripThumbURL({dimensions:400,region:reg,crs:'EPSG:3857'},(u,e)=>e?rej(e):res(u)));
  const r=await fetch(url); window.__blob=await r.blob(); window.__blobName='s2_'+name+'.png'; window.__save(); return r.status;
};
window.__res = window.__res || {};
