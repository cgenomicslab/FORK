// =====================================================================
// Shared utilities — used by all pages
// =====================================================================

function showEl(id) {
    const el = document.getElementById(id);
    if (el) el.classList.remove('hidden');
}

function hideEl(id) {
    const el = document.getElementById(id);
    if (el) el.classList.add('hidden');
}

function getEl(id) { return document.getElementById(id); }

function setHTML(id, html) {
    const el = getEl(id);
    if (el) el.innerHTML = html;
}

// Build an HTML table from an array of row objects
function makeTable(rows, columns) {
    if (!rows || rows.length === 0) return '<p class="text-muted">No data.</p>';
    const cols = columns || Object.keys(rows[0]);
    let html = '<div class="table-wrapper"><table><thead><tr>';
    cols.forEach(c => { html += `<th>${c}</th>`; });
    html += '</tr></thead><tbody>';
    rows.forEach(row => {
        html += '<tr>';
        cols.forEach(c => {
            const v = row[c];
            const display = v === null || v === undefined ? '' : String(v);
            if (c === 'accession' && display) {
                const url = 'https://www.uniprot.org/uniprotkb/' + encodeURIComponent(display) + '/entry';
                html += `<td title="View ${display} on UniProt"><a href="${url}" target="_blank" rel="noopener">${display}</a></td>`;
            } else {
                html += `<td title="${display}">${display}</td>`;
            }
        });
        html += '</tr>';
    });
    html += '</tbody></table></div>';
    return html;
}

// Render a base64 PNG in a container
function showImage(containerId, b64) {
    const el = getEl(containerId);
    if (!el) return;
    el.innerHTML = `<div class="result-img-wrap"><img src="data:image/png;base64,${b64}"></div>`;
    el.classList.remove('hidden');
}

// Show an alert message
function showAlert(containerId, type, message) {
    const el = getEl(containerId);
    if (!el) return;
    el.innerHTML = `<div class="alert alert-${type}">${message}</div>`;
    el.classList.remove('hidden');
}

function clearAlert(containerId) {
    const el = getEl(containerId);
    if (el) { el.innerHTML = ''; el.classList.add('hidden'); }
}

// Trigger a file download from a string
function downloadText(filename, content, mime) {
    const blob = new Blob([content], { type: mime || 'text/plain' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    a.click();
    URL.revokeObjectURL(a.href);
}

// Parse comma/newline separated IDs into an array of integers
function parseTaxIds(str) {
    if (!str || !str.trim()) return null;
    const ids = [];
    str.split(/[\s,]+/).forEach(tok => {
        const n = parseInt(tok.trim(), 10);
        if (!isNaN(n)) ids.push(n);
    });
    return ids.length ? ids : null;
}

// Poll a job until done
function pollJob(url, interval, onProgress, onDone, onError) {
    const timer = setInterval(async () => {
        try {
            const res  = await fetch(url);
            const data = await res.json();
            if (data.log && data.log.length) onProgress(data.log);
            if (data.status === 'done') {
                clearInterval(timer);
                onDone(data);
            } else if (data.status === 'error') {
                clearInterval(timer);
                onError(data.error || 'Unknown error');
            }
        } catch (e) {
            clearInterval(timer);
            onError(e.message);
        }
    }, interval);
    return timer;
}

// Persist job IDs across page navigations
function saveJobId(key, jobId) { localStorage.setItem('job_' + key, jobId); }
function loadJobId(key)        { return localStorage.getItem('job_' + key); }
function clearJobId(key)       { localStorage.removeItem('job_' + key); }

// Tab switcher
function initTabs(containerSelector) {
    const container = document.querySelector(containerSelector);
    if (!container) return;
    const buttons  = container.querySelectorAll('.tab-btn');
    const contents = container.querySelectorAll('.tab-content');
    buttons.forEach(btn => {
        btn.addEventListener('click', () => {
            buttons.forEach(b  => b.classList.remove('active'));
            contents.forEach(c => c.classList.remove('active'));
            btn.classList.add('active');
            const target = container.querySelector('#' + btn.dataset.tab);
            if (target) target.classList.add('active');
        });
    });
    if (buttons.length) buttons[0].click();
}

// =====================================================================
// =====================================================================
// Custom file-upload button  (replaces "Περιήγηση" / "Browse")
// =====================================================================

function initFileInputs() {
    document.querySelectorAll('input[type="file"].file-input').forEach(input => {
        if (input.dataset.customized) return;
        input.dataset.customized = '1';

        // Build the replacement row
        const row = document.createElement('div');
        row.className = 'file-upload-row';
        input.parentNode.insertBefore(row, input);

        const label = document.createElement('label');
        label.setAttribute('for', input.id);
        label.className = 'btn btn-secondary btn-sm file-upload-label';
        label.textContent = 'Upload';

        const display = document.createElement('span');
        display.className = 'file-name-display';
        display.textContent = 'No file selected';

        // Move the input into the row and hide the browser UI
        input.classList.add('file-input-hidden');
        row.appendChild(label);
        row.appendChild(display);
        row.appendChild(input);

        input.addEventListener('change', () => {
            display.textContent = input.files[0]?.name || 'No file selected';
        });
    });
}

document.addEventListener('DOMContentLoaded', initFileInputs);

// =====================================================================
// Active nav link
// =====================================================================

document.addEventListener('DOMContentLoaded', () => {
    const path = window.location.pathname;
    document.querySelectorAll('.main-nav a, .sidebar a.nav-item').forEach(a => {
        if (a.getAttribute('href') === path) a.classList.add('active');
    });
});

// =====================================================================
// DB Config Panel
// =====================================================================

let _dbPanelOpen = false;

function toggleDbPanel() {
    _dbPanelOpen = !_dbPanelOpen;
    const panel = getEl('db-panel');
    if (panel) panel.classList.toggle('open', _dbPanelOpen);
}

function closeDbPanel() {
    _dbPanelOpen = false;
    const panel = getEl('db-panel');
    if (panel) panel.classList.remove('open');
}

async function applyDbConfig() {
    const host     = (getEl('db-host')     || {}).value || '';
    const user     = (getEl('db-user')     || {}).value || '';
    const password = (getEl('db-password') || {}).value || '';
    const database = (getEl('db-name')     || {}).value || '';
    const port     = (getEl('db-port')     || {}).value || '';

    const statusEl = getEl('db-panel-status');
    if (statusEl) statusEl.textContent = 'Saving…';

    try {
        await fetch('/api/db-config', {
            method:  'POST',
            headers: { 'Content-Type': 'application/json' },
            body:    JSON.stringify({ host, user, password, database, port }),
        });
        if (statusEl) statusEl.textContent = 'Saved. Testing connection…';
        await checkDbStatus();
        closeDbPanel();
    } catch (e) {
        if (statusEl) statusEl.textContent = 'Error: ' + e.message;
    }
}

// Check DB connectivity and update the header indicator
async function checkDbStatus() {
    const dot  = getEl('header-db-dot');
    if (dot) dot.className = 'db-dot checking';

    try {
        const res  = await fetch('/api/db-info');
        const data = await res.json();
        if (data.error) {
            if (dot) dot.className = 'db-dot error';
            return false;
        }
        if (dot) dot.className = 'db-dot connected';
        return true;
    } catch (e) {
        if (dot) dot.className = 'db-dot error';
        return false;
    }
}

// Pre-populate form from .env defaults, then check connection
async function initDbConfig() {
    try {
        const res  = await fetch('/api/db-defaults');
        const data = await res.json();

        const hostEl = getEl('db-host');
        const userEl = getEl('db-user');
        const nameEl = getEl('db-name');

        // Only pre-fill if the fields are empty
        const portEl = getEl('db-port');
        if (hostEl && !hostEl.value) hostEl.value = data.host     || '';
        if (userEl && !userEl.value) userEl.value = data.user     || '';
        if (nameEl && !nameEl.value) nameEl.value = data.database || '';
        if (portEl && !portEl.value) portEl.value = data.port     || '3306';
        // Never pre-fill the password field for security
    } catch (_) { /* ignore — env might not have defaults */ }

    checkDbStatus();
}

// Close DB panel when clicking outside
document.addEventListener('click', e => {
    if (!_dbPanelOpen) return;
    const panel  = getEl('db-panel');
    const toggle = getEl('db-toggle-btn');
    if (panel && !panel.contains(e.target) && toggle && !toggle.contains(e.target)) {
        closeDbPanel();
    }
});

// Wire up DB panel buttons
document.addEventListener('DOMContentLoaded', () => {
    const toggleBtn = getEl('db-toggle-btn');
    const cancelBtn = getEl('db-cancel-btn');
    const applyBtn  = getEl('db-apply-btn');

    if (toggleBtn) toggleBtn.addEventListener('click', toggleDbPanel);
    if (cancelBtn) cancelBtn.addEventListener('click', closeDbPanel);
    if (applyBtn)  applyBtn.addEventListener('click', applyDbConfig);

    initDbConfig();
});
