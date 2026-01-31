// Auto-detect language direction based on URL path
// Arabic is primary (root level) - RTL
// English is secondary (/en/) - LTR

(function() {
  'use strict';

  // Run on page load
  function setDirection() {
    const path = window.location.pathname;
    const html = document.documentElement;

    // Check if we're on an English page (/en/)
    if (path.includes('/en/') || path.startsWith('/en')) {
      // English pages - LTR
      html.setAttribute('dir', 'ltr');
      html.setAttribute('lang', 'en');
      document.body.setAttribute('dir', 'ltr');
    } else {
      // Default to RTL for Arabic (root level, /blog/, /courses/)
      html.setAttribute('dir', 'rtl');
      html.setAttribute('lang', 'ar');
      document.body.setAttribute('dir', 'rtl');
    }
  }

  // Set direction immediately
  setDirection();

  // Also set on navigation (for SPA-like navigation in Material)
  if (document.addEventListener) {
    document.addEventListener('DOMContentLoaded', setDirection);

    // Material for MkDocs uses instant navigation
    document.addEventListener('DOMContentLoaded', function() {
      // Re-run on any navigation event
      const observer = new MutationObserver(function() {
        setDirection();
      });

      observer.observe(document.body, {
        childList: true,
        subtree: true
      });
    });
  }
})();
