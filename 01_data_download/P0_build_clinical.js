// P0 v2: build definitive TARGET-OS clinical table from expanded GDC cases JSON
// Outputs: data/processed/TARGET-OS_clinical.csv + summary stats
const fs = require('fs');
const root = 'D:/projects/OS_matrix_dualKO';
const j = JSON.parse(fs.readFileSync(root + '/data/raw/TARGET-OS/gdc_cases_expanded.json', 'utf8'));

const EFS_EVENTS = new Set(['Relapse', 'Progression', 'Death', 'Death without Remission',
  'Second Malignant Neoplasm', 'Presented with Metastases']);
const MET_EVENT = 'Presented with Metastases';

const rows = [];
j.data.hits.forEach(h => {
  const d = h.demographic || {};
  const dx = (h.diagnoses || [])[0] || {};
  const fus = h.follow_ups || [];
  const vital = d.vital_status || '';

  // ---- OS ----
  const lastContacts = fus.filter(f => f.timepoint_category === 'Last Contact')
    .map(f => f.days_to_follow_up).filter(v => v != null);
  const dxLastFu = dx.days_to_last_follow_up;
  let osTime = null;
  if (vital === 'Dead') osTime = d.days_to_death;
  else osTime = Math.max(dxLastFu != null ? dxLastFu : -1, ...(lastContacts.length ? lastContacts : [-1]));
  if (osTime === -1 || osTime == null) osTime = '';
  const osEvent = vital === 'Dead' ? 1 : (vital === 'Alive' ? 0 : '');

  // ---- EFS ----
  const evFus = fus.filter(f => f.first_event != null && EFS_EVENTS.has(f.first_event)
    && (f.days_to_first_event != null || f.days_to_follow_up != null));
  let efsTime = '', efsEvent = '';
  if (vital === 'Dead' && d.days_to_death != null && !evFus.length) {
    efsTime = d.days_to_death; efsEvent = 1;  // death without recorded relapse
  } else if (evFus.length) {
    const first = evFus.reduce((a, b) => ((a.days_to_first_event ?? a.days_to_follow_up) <= (b.days_to_first_event ?? b.days_to_follow_up) ? a : b));
    efsTime = first.days_to_first_event ?? first.days_to_follow_up;
    efsEvent = 1;
  } else {
    const cens = Math.max(dxLastFu != null ? dxLastFu : -1, ...(lastContacts.length ? lastContacts : [-1]));
    if (cens > 0) { efsTime = cens; efsEvent = 0; }
  }

  // ---- metastasis at presentation ----
  const metPresent = fus.some(f => f.first_event === MET_EVENT) ? 1 : 0;

  rows.push({
    submitter_id: h.submitter_id,
    gender: d.sex_at_birth || d.gender || '',
    vital_status: vital,
    age_at_diagnosis_years: dx.age_at_diagnosis != null ? (dx.age_at_diagnosis / 365.25).toFixed(2) : '',
    primary_site: dx.primary_diagnosis || h.primary_site || '',
    metastasis_at_presentation: metPresent,
    days_to_last_follow_up: dxLastFu ?? '',
    OS_time_days: osTime, OS_event: osEvent,
    EFS_time_days: efsTime, EFS_event: efsEvent,
    first_event_types: [...new Set(fus.map(f => f.first_event).filter(Boolean))].join('|')
  });
});
rows.sort((a, b) => a.submitter_id.localeCompare(b.submitter_id));
const cols = Object.keys(rows[0]);
const csv = [cols.join(',')].concat(rows.map(r => cols.map(c => String(r[c]).replace(/,/g, ';')).join(','))).join('\n');
fs.writeFileSync(root + '/data/processed/TARGET-OS_clinical.csv', csv);

const os = rows.filter(r => r.OS_time_days !== '' && r.OS_event !== '');
const efs = rows.filter(r => r.EFS_time_days !== '' && r.EFS_event !== '');
console.log('total cases:', rows.length);
console.log('OS usable:', os.length, '(deaths:', os.filter(r => r.OS_event === 1).length, ', censored:', os.filter(r => r.OS_event === 0).length + ')');
console.log('EFS usable:', efs.length, '(events:', efs.filter(r => r.EFS_event === 1).length, ', censored:', efs.filter(r => r.EFS_event === 0).length + ')');
console.log('metastasis at presentation:', rows.filter(r => r.metastasis_at_presentation === 1).length);
