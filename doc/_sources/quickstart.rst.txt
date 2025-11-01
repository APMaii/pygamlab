Getting Started
===============

.. raw:: html

    <style>
    /* ======== Global Layout ======== */
    .getting-started-wrapper {
        display: flex;
        justify-content: center;
        align-items: center;
        flex-direction: column;
        margin: 3rem auto 4rem auto;
        max-width: 800px;
        padding: 2rem;
    }

    /* ======== Glassmorphic Card ======== */
    .getting-started-card {
        position: relative;
        border-radius: 18px;
        padding: 2.5rem 2rem 2.8rem 2rem;
        text-align: center;
        backdrop-filter: blur(14px);
        -webkit-backdrop-filter: blur(14px);
        box-shadow: 0 8px 24px rgba(0, 0, 0, 0.15);
        border: 1px solid rgba(255, 255, 255, 0.15);
        transition: transform 0.3s ease, box-shadow 0.3s ease;
    }

    .getting-started-card:hover {
        transform: translateY(-6px);
        box-shadow: 0 12px 28px rgba(255, 106, 0, 0.25);
    }

    /* ======== Color Themes ======== */
    html[data-theme="light"] .getting-started-card {
        background: linear-gradient(160deg, #fff9f4 0%, #fff3e6 100%);
        color: #1b1b1b;
    }

    html[data-theme="dark"] .getting-started-card {
        background: linear-gradient(160deg, #1a120c 0%, #2b1f14 100%);
        color: #f4f4f4;
        border: 1px solid rgba(255, 106, 0, 0.25);
    }

    /* ======== Header ======== */
    .getting-started-card h2 {
        font-size: 2rem;
        margin-bottom: 0.5rem;
        color: #ff6a00;
        font-weight: 800;
        letter-spacing: 0.3px;
    }

    .getting-started-card p {
        font-size: 1.05rem;
        line-height: 1.7;
        margin-bottom: 1.3rem;
    }

    /* ======== Feature List ======== */
    .getting-started-list {
        text-align: left;
        display: inline-block;
        margin: 1.6rem 0;
        font-size: 0.98rem;
        line-height: 1.6;
    }

    .getting-started-list li {
        margin-bottom: 0.4rem;
        list-style: none;
        position: relative;
        padding-left: 1.6rem;
    }

    .getting-started-list li::before {
        content: "🔥";
        position: absolute;
        left: 0;
        top: 0.05rem;
    }

    /* ======== CTA Button ======== */
    .getting-started-btn {
        display: inline-block;
        margin-top: 1.4rem;
        background: linear-gradient(90deg, #ff6a00, #ff8f38);
        color: #fff;
        padding: 0.85rem 2rem;
        font-size: 1.05rem;
        border-radius: 10px;
        text-decoration: none;
        font-weight: 600;
        letter-spacing: 0.3px;
        box-shadow: 0 4px 12px rgba(255,106,0,0.4);
        transition: all 0.3s ease;
    }

    .getting-started-btn:hover {
        background: linear-gradient(90deg, #ff7f1a, #ffaa52);
        transform: translateY(-4px);
        box-shadow: 0 8px 20px rgba(255,106,0,0.55);
    }

    /* ======== Subtle Glow Ring ======== */
    .glow-ring {
        position: absolute;
        top: -1px;
        left: -1px;
        right: -1px;
        bottom: -1px;
        border-radius: 18px;
        background: linear-gradient(120deg, rgba(255,106,0,0.25), rgba(255,255,255,0.1));
        z-index: -1;
    }
    </style>

    <div class="getting-started-wrapper">
      <div class="getting-started-card">
        <div class="glow-ring"></div>
        <h2>🚀 Getting Started with PyGamLab</h2>
        <p><b>PyGamLab</b> is a scientific Python library for researchers, engineers, and students in
        <b>materials science</b> and <b>nanotechnology</b>, providing tools for structure design,
        AI integration, and data-driven discovery.</p>

        <p>Begin your journey by exploring the comprehensive Quickstart tutorial:</p>

        <ul class="getting-started-list">
          <li>⚙️ Installation instructions</li>
          <li>📦 Core modules overview</li>
          <li>📈 Step-by-step examples</li>
          <li>🧠 API reference</li>
        </ul>

        <a href="quickstart_tutorial.html" class="getting-started-btn">
          Start Quickstart Tutorial →
        </a>
      </div>
    </div>
