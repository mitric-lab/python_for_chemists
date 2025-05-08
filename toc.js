// Populate the sidebar
//
// This is a script, and not included directly in the page, to control the total size of the book.
// The TOC contains an entry for each page, so if each page includes a copy of the TOC,
// the total size of the page becomes O(n**2).
class MDBookSidebarScrollbox extends HTMLElement {
    constructor() {
        super();
    }
    connectedCallback() {
        this.innerHTML = '<ol class="chapter"><li class="chapter-item expanded affix "><li class="part-title">Inhaltsverzeichnis</li><li class="chapter-item expanded "><a href="00-preface.html"><strong aria-hidden="true">0.</strong> Vorwort</a><a class="toggle"><div>❱</div></a></li><li><ol class="section"><li class="chapter-item "><a href="00-preface/01-motivation.html"><strong aria-hidden="true">0.1.</strong> Motivation</a></li><li class="chapter-item "><a href="00-preface/02-getting_started.html"><strong aria-hidden="true">0.2.</strong> Erste Schritte</a></li><li class="chapter-item "><a href="00-preface/03-mdbook_usage.html"><strong aria-hidden="true">0.3.</strong> Bedienung dieser Website</a></li></ol></li><li class="chapter-item expanded "><li class="spacer"></li><li class="chapter-item expanded "><a href="01-regression.html"><strong aria-hidden="true">1.</strong> Regression Analysis</a><a class="toggle"><div>❱</div></a></li><li><ol class="section"><li class="chapter-item "><a href="01-regression/01-least_squares.html"><strong aria-hidden="true">1.1.</strong> Least Squares</a></li><li class="chapter-item "><a href="01-regression/02-linear_regression.html"><strong aria-hidden="true">1.2.</strong> Linear Regression</a></li><li class="chapter-item "><a href="01-regression/03-numerical_optimisation.html"><strong aria-hidden="true">1.3.</strong> Numerical Optimisation</a></li><li class="chapter-item "><div><strong aria-hidden="true">1.4.</strong> Nonlinear Regression</div></li></ol></li><li class="chapter-item expanded "><div><strong aria-hidden="true">2.</strong> Differentialgleichungen</div><a class="toggle"><div>❱</div></a></li><li><ol class="section"><li class="chapter-item "><div><strong aria-hidden="true">2.1.</strong> Anfangswertproblem</div></li><li class="chapter-item "><div><strong aria-hidden="true">2.2.</strong> Euler-Verfahren</div></li><li class="chapter-item "><div><strong aria-hidden="true">2.3.</strong> Runge-Kutta-Verfahren</div></li><li class="chapter-item "><div><strong aria-hidden="true">2.4.</strong> Finite-Differenzen-Verfahren</div></li></ol></li><li class="chapter-item expanded "><div><strong aria-hidden="true">3.</strong> Eigenwert- und Singulärwertzerlegung</div><a class="toggle"><div>❱</div></a></li><li><ol class="section"><li class="chapter-item "><div><strong aria-hidden="true">3.1.</strong> Eigenwertzerlegung</div></li><li class="chapter-item "><div><strong aria-hidden="true">3.2.</strong> Singulärwertzerlegung</div></li><li class="chapter-item "><div><strong aria-hidden="true">3.3.</strong> Hauptkomponentenanalyse</div></li><li class="chapter-item "><div><strong aria-hidden="true">3.4.</strong> Hauptkoordinatenanalyse</div></li><li class="chapter-item "><div><strong aria-hidden="true">3.5.</strong> Lineare Gleichungssysteme</div></li></ol></li><li class="chapter-item expanded "><div><strong aria-hidden="true">4.</strong> Maschinelles Lernen</div><a class="toggle"><div>❱</div></a></li><li><ol class="section"><li class="chapter-item "><div><strong aria-hidden="true">4.1.</strong> Überwachtes Lernen</div></li><li class="chapter-item "><div><strong aria-hidden="true">4.2.</strong> Unüberwachtes Lernen</div></li></ol></li><li class="chapter-item expanded "><div><strong aria-hidden="true">5.</strong> Neuronale Netzwerke</div><a class="toggle"><div>❱</div></a></li><li><ol class="section"><li class="chapter-item "><div><strong aria-hidden="true">5.1.</strong> Single-Layer-Perzeptron</div></li><li class="chapter-item "><div><strong aria-hidden="true">5.2.</strong> Multi-Layer-Perzeptron</div></li></ol></li><li class="chapter-item expanded "><li class="spacer"></li><li class="chapter-item expanded affix "><div>Zusammenfassung und Ausblick</div></li><li class="chapter-item expanded affix "><li class="spacer"></li><li class="chapter-item expanded affix "><a href="psets/01.html">Problem Set 1</a></li><li class="chapter-item expanded affix "><div>Problem Set 2</div></li><li class="chapter-item expanded affix "><div>Problem Set 3</div></li><li class="chapter-item expanded affix "><div>Problem Set 4</div></li><li class="chapter-item expanded affix "><div>Problem Set 5</div></li><li class="chapter-item expanded affix "><div>Problem Set 6</div></li><li class="chapter-item expanded affix "><div>Sample Exam</div></li></ol>';
        // Set the current, active page, and reveal it if it's hidden
        let current_page = document.location.href.toString().split("#")[0];
        if (current_page.endsWith("/")) {
            current_page += "index.html";
        }
        var links = Array.prototype.slice.call(this.querySelectorAll("a"));
        var l = links.length;
        for (var i = 0; i < l; ++i) {
            var link = links[i];
            var href = link.getAttribute("href");
            if (href && !href.startsWith("#") && !/^(?:[a-z+]+:)?\/\//.test(href)) {
                link.href = path_to_root + href;
            }
            // The "index" page is supposed to alias the first chapter in the book.
            if (link.href === current_page || (i === 0 && path_to_root === "" && current_page.endsWith("/index.html"))) {
                link.classList.add("active");
                var parent = link.parentElement;
                if (parent && parent.classList.contains("chapter-item")) {
                    parent.classList.add("expanded");
                }
                while (parent) {
                    if (parent.tagName === "LI" && parent.previousElementSibling) {
                        if (parent.previousElementSibling.classList.contains("chapter-item")) {
                            parent.previousElementSibling.classList.add("expanded");
                        }
                    }
                    parent = parent.parentElement;
                }
            }
        }
        // Track and set sidebar scroll position
        this.addEventListener('click', function(e) {
            if (e.target.tagName === 'A') {
                sessionStorage.setItem('sidebar-scroll', this.scrollTop);
            }
        }, { passive: true });
        var sidebarScrollTop = sessionStorage.getItem('sidebar-scroll');
        sessionStorage.removeItem('sidebar-scroll');
        if (sidebarScrollTop) {
            // preserve sidebar scroll position when navigating via links within sidebar
            this.scrollTop = sidebarScrollTop;
        } else {
            // scroll sidebar to current active section when navigating via "next/previous chapter" buttons
            var activeSection = document.querySelector('#sidebar .active');
            if (activeSection) {
                activeSection.scrollIntoView({ block: 'center' });
            }
        }
        // Toggle buttons
        var sidebarAnchorToggles = document.querySelectorAll('#sidebar a.toggle');
        function toggleSection(ev) {
            ev.currentTarget.parentElement.classList.toggle('expanded');
        }
        Array.from(sidebarAnchorToggles).forEach(function (el) {
            el.addEventListener('click', toggleSection);
        });
    }
}
window.customElements.define("mdbook-sidebar-scrollbox", MDBookSidebarScrollbox);
