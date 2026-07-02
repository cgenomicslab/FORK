# Hosting Guide — UniProt Lab Manager

This guide explains how to host the Flask app on your own server so other people in your lab (or outside) can use it from a browser without installing anything.

---

## How a web app works (the basics)

When you run `python app.py` on your laptop, only you can access it at `localhost:5000`. To let others access it, you need to:

1. Run it on a machine that is always on (your server)
2. Make sure that machine is reachable over the network (open port, domain name)
3. Keep the app running even when no one is using it (process manager)
4. Optionally: put something in front of it (a reverse proxy) to handle HTTPS, multiple users, and clean URLs

---

## The four options

---

## Option 1 — Flask development server (simplest, not recommended for others)

### What it is
Just run the app directly with Python, the same way you do on your laptop.

### How to run it
```bash
# on the server
python app.py
```
Or, to make it accessible from outside the server:
```bash
flask run --host=0.0.0.0 --port=5000
```
Users go to `http://your-server-ip:5000`

### Infrastructure
```
User browser  →  your-server-ip:5000  →  Flask (Python)
```

### Pros
- Zero setup
- Good for testing that everything works on the server

### Cons
- Flask's built-in server handles **one request at a time** — if two people click at the same moment, one waits
- Crashes if there's an error, and doesn't restart itself
- No HTTPS
- Not meant for real use — Flask itself warns you "do not use in production"
- If you close your terminal / SSH session, the app stops

### When to use it
Only for quick testing on the server. Never as a real deployment.

---

## Option 2 — Gunicorn (production Python server)

### What it is
Gunicorn is a proper Python web server. It replaces Flask's built-in server and can handle many users at the same time by running multiple worker processes.

### How to install
```bash
pip install gunicorn
# or via conda:
conda install -c conda-forge gunicorn
```

### How to run it
```bash
gunicorn -w 4 -b 0.0.0.0:5000 app:app
```
- `-w 4` → 4 worker processes (can handle 4 requests simultaneously; set to `2 × CPU cores + 1`)
- `-b 0.0.0.0:5000` → listen on all network interfaces, port 5000
- `app:app` → the file is `app.py`, the Flask object inside it is also called `app`

### Infrastructure
```
User browser  →  your-server-ip:5000  →  Gunicorn (4 workers)  →  your Flask app
```

### Pros
- Handles multiple users at the same time
- Stable and battle-tested
- Easy to install (one pip command)
- Industry standard for Flask/Django apps

### Cons
- Still no HTTPS on its own
- Still need to manage the process (it won't restart on reboot unless you add a process manager — see Option 4)
- Port 5000 is exposed directly, which is fine for a private lab server but not ideal for public access

### When to use it
Good starting point for lab-internal use. Combine with Option 4 (Systemd) to keep it running permanently.

---

## Option 3 — Gunicorn + Nginx (recommended for real deployment)

### What it is
Nginx is a lightweight web server that sits **in front** of Gunicorn. It handles the "public-facing" side: HTTPS, domain names, serving static files, and forwarding requests to Gunicorn. This is called a **reverse proxy**.

### Why add Nginx?
- Gunicorn is good at running Python, but slow at serving static files (CSS, images, JS). Nginx is extremely fast at this.
- Nginx handles SSL/HTTPS certificates
- You can host multiple apps on the same server on port 80/443 (one Nginx, multiple backends)
- Adds a security layer between the internet and your Python process

### How the request flows
```
User browser
     │
     ▼  port 443 (HTTPS) or 80 (HTTP)
   Nginx
     │  proxies the request internally
     ▼  port 5000 (only reachable on the server itself)
  Gunicorn
     │
     ▼
  Flask app
```

### Nginx config (basic example)
Create a file at `/etc/nginx/sites-available/uniprot-lab-manager`:

```nginx
server {
    listen 80;
    server_name yourdomain.com;   # or your server IP

    # Serve static files directly (fast, no Python involved)
    location /static/ {
        alias /path/to/uniprot-lab-manager-copy/static/;
    }

    # Forward everything else to Gunicorn
    location / {
        proxy_pass http://127.0.0.1:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    }
}
```

Enable it:
```bash
sudo ln -s /etc/nginx/sites-available/uniprot-lab-manager /etc/nginx/sites-enabled/
sudo nginx -t        # test config
sudo systemctl reload nginx
```

### Adding HTTPS (free, with Let's Encrypt)
If your server has a public domain name:
```bash
sudo apt install certbot python3-certbot-nginx
sudo certbot --nginx -d yourdomain.com
```
Certbot automatically edits your Nginx config and renews the certificate every 90 days.

### Pros
- Proper HTTPS with a real certificate
- Fast static file serving
- Can host multiple apps on one server
- Industry standard — what most labs and companies use
- Nginx is extremely stable and low-resource

### Cons
- More to set up (two things to configure instead of one)
- Need a domain name for HTTPS (if you only have an IP, you can use a self-signed certificate but browsers will warn users)

### When to use it
Whenever you want a clean, stable deployment that others outside the server can reach. Especially if you have a domain name.

---

## Option 4 — Systemd (keep the app alive permanently)

### What it is
Systemd is Linux's built-in service manager. It can start your app automatically when the server boots, and restart it automatically if it crashes. This is **not an alternative** to the above options — it works alongside them.

### Create a service file
Create `/etc/systemd/system/uniprot-lab-manager.service`:

```ini
[Unit]
Description=UniProt Lab Manager Flask App
After=network.target

[Service]
User=your_linux_username
WorkingDirectory=/path/to/uniprot-lab-manager-copy
EnvironmentFile=/path/to/uniprot-lab-manager-copy/.env
ExecStart=/path/to/conda/envs/bio_tools/bin/gunicorn -w 4 -b 127.0.0.1:5000 app:app
Restart=always
RestartSec=5

[Install]
WantedBy=multi-user.target
```

Enable and start:
```bash
sudo systemctl daemon-reload
sudo systemctl enable uniprot-lab-manager   # start on boot
sudo systemctl start uniprot-lab-manager    # start now
sudo systemctl status uniprot-lab-manager   # check it's running
```

Useful commands:
```bash
sudo systemctl stop uniprot-lab-manager     # stop
sudo systemctl restart uniprot-lab-manager  # restart (e.g. after code update)
journalctl -u uniprot-lab-manager -f        # live logs
```

### Pros
- App survives reboots and crashes automatically
- Clean logs via `journalctl`
- Standard Linux way — no extra software needed

### Cons
- Requires root/sudo on the server to set up
- Linux-only (but your server is almost certainly Linux)

---

## Option 5 — Docker

### What it is
Docker packages the entire app — Python, dependencies, config — into a self-contained **container**. The container runs identically on any machine that has Docker installed, regardless of what else is on the server.

### How it works
You write a `Dockerfile` that describes how to build the environment:

```dockerfile
FROM continuumio/miniconda3

WORKDIR /app
COPY uniprot-lab-manager.yml .
RUN conda env create -f uniprot-lab-manager.yml

COPY . .

EXPOSE 5000
CMD ["conda", "run", "-n", "bio_tools", "gunicorn", "-w", "4", "-b", "0.0.0.0:5000", "app:app"]
```

Build and run:
```bash
docker build -t uniprot-lab-manager .
docker run -d -p 5000:5000 --env-file .env uniprot-lab-manager
```

Or with Docker Compose (easier to manage with a database too):
```yaml
# docker-compose.yml
services:
  app:
    build: .
    ports:
      - "5000:5000"
    env_file: .env
    restart: always
```
```bash
docker compose up -d
```

### Pros
- Completely isolated from the rest of the server — no dependency conflicts
- Easy to move to a different server (just copy the image)
- Easy to update: rebuild image, restart container
- Can version your deployments
- Docker Compose makes it easy to also run the MySQL database in a container alongside the app

### Cons
- Docker itself needs to be installed and learned
- Slightly more overhead than running natively
- With Conda inside Docker, images can be large (1–2 GB)
- Your MySQL database needs to be accessible from inside the container (network config)

### When to use it
When you want clean, reproducible deployments or plan to move the app between servers. Also good if you want to containerize the MySQL database too.

---

## Comparison table

| | Option 1 (Flask dev) | Option 2 (Gunicorn) | Option 3 (+ Nginx) | Option 4 (+ Systemd) | Option 5 (Docker) |
|---|---|---|---|---|---|
| Multiple users | ❌ | ✅ | ✅ | ✅ | ✅ |
| HTTPS | ❌ | ❌ | ✅ | ✅ (with Nginx) | ✅ (with Nginx) |
| Survives reboot | ❌ | ❌ | ❌ | ✅ | ✅ |
| Auto-restart on crash | ❌ | ❌ | ❌ | ✅ | ✅ |
| Setup difficulty | Very easy | Easy | Medium | Medium | Medium–Hard |
| Good for lab use | Testing only | Internal only | ✅ Yes | ✅ Yes | ✅ Yes |

---

## What I recommend for your lab

**Short answer: Gunicorn + Nginx + Systemd**

It is the standard stack for exactly this kind of lab tool. Here is the setup in plain terms:

- **Gunicorn** runs your Flask app and handles multiple users
- **Nginx** sits in front, serves your CSS/images fast, and handles HTTPS if you have a domain
- **Systemd** keeps everything running 24/7 and restarts after reboots

The full setup takes about 30–60 minutes the first time.

**If you also want the MySQL database to be easy to manage alongside the app**, consider adding Docker Compose later — it lets you run both the app and the database as containers with one command.

---

## Things to check on your server before starting

1. **What OS?** — These instructions assume Linux (Ubuntu/Debian). Let me know if it's different.
2. **Is port 80/443 open?** — Check with your IT/sysadmin if you're behind a firewall.
3. **Do you have a domain name?** — Needed for proper HTTPS. If not, you can use the server's IP with a self-signed certificate (browsers will show a warning).
4. **Where is the MySQL database?** — If it's on the same server, the `.env` file stays as-is. If it's on a different machine, update `DB_HOST`.
5. **Conda environment** — Make sure the `bio_tools` conda environment is installed on the server and `gunicorn` is added to it.

---

## Updating the app after code changes

With Gunicorn + Systemd, deploying an update is:

```bash
cd /path/to/uniprot-lab-manager-copy
git pull
sudo systemctl restart uniprot-lab-manager
```

That's it — three commands.
