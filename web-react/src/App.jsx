import { lazy, Suspense, useEffect, useRef, useState } from 'react';
import {
  AlertCircle,
  AlertTriangle,
  Aperture,
  Beaker,
  BookOpen,
  ChevronDown,
  Database,
  Download,
  FileText,
  GitCommit,
  Info,
  ListChecks,
  ShieldCheck,
  Sparkles,
} from 'lucide-react';

const SpectralHero = lazy(() => import('./SpectralHero.jsx'));

const warningCategories = [
  {
    key: 'degenerate-split',
    title: 'Degenerate comparison split',
    description: 'The released flag values place every selected object in one cohort, so a two-group comparison is not defined.',
    test: (warning) => warning.includes('split is degenerate'),
    tone: 'caveat',
  },
  {
    key: 'empty-flagged-data',
    title: 'Flagged-cohort measurements unavailable',
    description: 'Group-specific statistics were skipped because the flagged cohort contains no finite measurements.',
    test: (warning) => warning.includes("group 'flagged'") || warning.includes("for group 'flagged'"),
    tone: 'quality',
  },
  {
    key: 'negative-control',
    title: 'Negative control unavailable',
    description: 'The permutation test requires at least two members in each cohort and was therefore not evaluated.',
    test: (warning) => warning.includes('negative control skipped'),
    tone: 'caveat',
  },
];

function parseResultJson(text) {
  const validJson = text.replace(/:\s*NaN(?=\s*[,}])/g, ': null');
  return JSON.parse(validJson);
}

function useJson(path) {
  const [state, setState] = useState({ data: null, error: null, loading: true });
  useEffect(() => {
    let cancelled = false;
    fetch(path)
      .then((response) => {
        if (!response.ok) throw new Error(`${path}: HTTP ${response.status}`);
        return response.text();
      })
      .then(parseResultJson)
      .then((data) => {
        if (!cancelled) setState({ data, error: null, loading: false });
      })
      .catch((error) => {
        if (!cancelled) setState({ data: null, error, loading: false });
      });
    return () => { cancelled = true; };
  }, [path]);
  return state;
}

function useNearViewport() {
  const ref = useRef(null);
  const [visible, setVisible] = useState(false);

  useEffect(() => {
    if (!ref.current || !('IntersectionObserver' in window)) {
      setVisible(true);
      return undefined;
    }
    const observer = new IntersectionObserver(
      ([entry]) => {
        if (entry.isIntersecting) {
          setVisible(true);
          observer.disconnect();
        }
      },
      { rootMargin: '180px' },
    );
    observer.observe(ref.current);
    return () => observer.disconnect();
  }, []);

  return [ref, visible];
}

function formatEstimate(estimate) {
  if (typeof estimate === 'number' && Number.isFinite(estimate)) return estimate.toPrecision(4);
  return 'Unavailable';
}

function Section({ icon: Icon, kicker, title, className = '', dark = false, children }) {
  return (
    <article className={`${dark ? 'spectral-panel-dark' : 'spectral-panel'} ${className}`}>
      <div className="mb-6 flex items-start gap-3">
        <span className={`section-mark ${dark ? 'section-mark-dark' : ''}`}>
          <Icon size={17} aria-hidden="true" />
        </span>
        <div>
          {kicker && <p className="spectral-kicker">{kicker}</p>}
          <h2 className={`font-display text-2xl leading-tight ${dark ? 'text-white' : 'text-plum-950'}`}>{title}</h2>
        </div>
      </div>
      {children}
    </article>
  );
}

function MetricCard({ metric, index }) {
  const hasUncertainty = metric.uncertainty_low != null && metric.uncertainty_high != null;
  const available = typeof metric.estimate === 'number' && Number.isFinite(metric.estimate);
  return (
    <article className={`metric-tile ${available ? '' : 'metric-tile-muted'}`}>
      <div className="flex items-center justify-between gap-3">
        <span className="font-mono text-[0.65rem] text-fuchsia-700">M{String(index + 1).padStart(2, '0')}</span>
        <span className="rounded-full border border-plum-200 bg-white/70 px-2 py-1 font-mono text-[0.62rem] text-plum-500">n={metric.sample_size}</span>
      </div>
      <p className="mt-5 break-words text-[0.68rem] uppercase leading-relaxed tracking-[0.12em] text-plum-600">
        {metric.name.replace(/_/g, ' ')}
      </p>
      <p className={`${available ? 'text-plum-950' : 'text-plum-400'} mt-3 font-display text-3xl leading-none`}>
        {formatEstimate(metric.estimate)}
      </p>
      <p className="mt-2 text-xs text-fuchsia-800/70">{metric.units}</p>
      {hasUncertainty && (
        <p className="mt-3 font-mono text-[0.67rem] text-plum-500">
          95% CI [{metric.uncertainty_low.toPrecision(3)}, {metric.uncertainty_high.toPrecision(3)}]
        </p>
      )}
    </article>
  );
}

function inverseNormalCDF(p) {
  if (p <= 0 || p >= 1) return NaN;
  const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
  const pLow = 0.02425;
  const pHigh = 1 - pLow;
  let q;
  let r;
  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
  if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  }
  q = Math.sqrt(-2 * Math.log(1 - p));
  return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
}

function ConfidenceExplorer({ metrics }) {
  const withCI = (metrics || []).filter((metric) => (
    typeof metric.estimate === 'number'
    && Number.isFinite(metric.estimate)
    && metric.uncertainty_low != null
    && metric.uncertainty_high != null
  ));
  const [selected, setSelected] = useState(null);
  const [confidence, setConfidence] = useState(95);

  useEffect(() => {
    if (!selected && withCI.length > 0) setSelected(withCI[0].name);
  }, [selected, withCI]);

  if (withCI.length === 0) return null;
  const metric = withCI.find((item) => item.name === selected) ?? withCI[0];
  const halfWidth95 = (metric.uncertainty_high - metric.uncertainty_low) / 2;
  const sigma = halfWidth95 / 1.959963984540054;
  const zLevel = inverseNormalCDF(0.5 + confidence / 200);
  const low = metric.estimate - zLevel * sigma;
  const high = metric.estimate + zLevel * sigma;

  return (
    <Section icon={Beaker} kicker="Interactive lens" title="Confidence-level explorer" className="h-full">
      <p className="text-sm leading-6 text-plum-700">
        Approximate interval derived from the reported 95% bootstrap bounds under a normal sampling
        distribution. The recorded 95% interval remains the computed result.
      </p>
      {withCI.length > 1 && (
        <select
          className="mt-5 w-full rounded-xl border border-plum-200 bg-white px-3 py-2 text-sm text-plum-900"
          value={metric.name}
          onChange={(event) => setSelected(event.target.value)}
        >
          {withCI.map((item) => <option key={item.name} value={item.name}>{item.name.replace(/_/g, ' ')}</option>)}
        </select>
      )}
      <label className="mt-6 flex items-center justify-between text-sm text-plum-700">
        <span>Confidence level</span>
        <span className="font-mono text-fuchsia-700">{confidence.toFixed(1)}%</span>
      </label>
      <input
        type="range"
        min="50"
        max="99.9"
        step="0.1"
        value={confidence}
        onChange={(event) => setConfidence(Number(event.target.value))}
        className="mt-3 w-full accent-fuchsia-700"
      />
      <p className="mt-5 font-display text-3xl text-plum-950">[{low.toPrecision(4)}, {high.toPrecision(4)}]</p>
      <p className="mt-2 text-xs text-plum-500">{metric.units} · estimate {metric.estimate.toPrecision(4)} · n={metric.sample_size}</p>
    </Section>
  );
}

function WarningAudit({ state }) {
  if (state.loading) return <p className="text-sm text-plum-500">Loading audit notes…</p>;
  if (state.error) {
    return (
      <div className="flex gap-3 rounded-xl border border-red-300 bg-red-50 p-4 text-sm text-red-900">
        <AlertCircle size={18} className="mt-0.5 shrink-0" aria-hidden="true" />
        Could not load results/warnings.json: {String(state.error)}
      </div>
    );
  }

  const entries = Array.isArray(state.data) ? state.data : [];
  if (entries.length === 0) {
    return (
      <div className="flex gap-3 rounded-xl border border-emerald-300 bg-emerald-50 p-4 text-sm text-emerald-900">
        <ShieldCheck size={18} className="mt-0.5 shrink-0" aria-hidden="true" />
        No warnings recorded in results/warnings.json.
      </div>
    );
  }

  const claimed = new Set();
  const groups = warningCategories.map((category) => {
    const items = entries.filter((warning, index) => {
      if (claimed.has(index) || !category.test(warning)) return false;
      claimed.add(index);
      return true;
    });
    return { ...category, items };
  }).filter((group) => group.items.length > 0);
  const unclassified = entries.filter((_, index) => !claimed.has(index));
  if (unclassified.length > 0) {
    groups.push({
      key: 'unclassified',
      title: 'Unclassified audit note',
      description: 'A record without a recognised presentation category; inspect the raw entry below.',
      tone: 'failure',
      items: unclassified,
    });
  }

  return (
    <div>
      <div className="mb-5 flex items-start gap-3 rounded-xl border border-fuchsia-200 bg-fuchsia-50 p-4">
        <Info size={18} className="mt-0.5 shrink-0 text-fuchsia-700" aria-hidden="true" />
        <p className="text-sm leading-6 text-plum-800">
          <strong>{entries.length} transparent audit notes</strong> explain why a flagged-cohort comparison
          is unavailable. This is the released null result, not a pipeline failure.
        </p>
      </div>
      <div className="grid gap-3 md:grid-cols-3">
        {groups.map((group) => (
          <div key={group.key} className={`warning-card warning-${group.tone}`}>
            <div className="flex items-start justify-between gap-3">
              <p className="font-display text-lg leading-tight text-plum-950">{group.title}</p>
              <span className="warning-count">{group.items.length}</span>
            </div>
            <p className="mt-3 text-xs leading-5 text-plum-600">{group.description}</p>
          </div>
        ))}
      </div>
      <details className="raw-warning-list mt-4 overflow-hidden rounded-xl border border-plum-200 bg-white">
        <summary className="flex cursor-pointer list-none items-center justify-between gap-3 p-4 text-sm font-semibold text-fuchsia-800">
          <span>Show all {entries.length} raw entries</span>
          <ChevronDown size={17} className="details-chevron" aria-hidden="true" />
        </summary>
        <ol className="space-y-3 border-t border-plum-100 p-5 pl-10 text-xs leading-5 text-plum-600">
          {entries.map((warning, index) => <li key={`${index}-${warning}`} className="pl-1 marker:text-fuchsia-500">{warning}</li>)}
        </ol>
      </details>
    </div>
  );
}

function LazySpectralHero() {
  const [ref, visible] = useNearViewport();
  return (
    <div ref={ref} className="spectral-hero-frame h-[22rem] overflow-hidden sm:h-[26rem] lg:h-full lg:min-h-[33rem]">
      {visible ? (
        <Suspense fallback={<div className="hero-loading grid h-full place-items-center text-xs uppercase tracking-[0.2em] text-fuchsia-200">Loading spectral instrument…</div>}>
          <SpectralHero />
        </Suspense>
      ) : (
        <div className="hero-loading h-full" aria-label="Spectral illustration placeholder" />
      )}
    </div>
  );
}

function CohortStatement({ metrics }) {
  const byName = Object.fromEntries(metrics.map((metric) => [metric.name, metric]));
  const clean = byName.n_clean?.estimate;
  const flagged = byName.n_flagged?.estimate;
  if (!Number.isFinite(clean) || !Number.isFinite(flagged)) return null;
  return (
    <div className="cohort-statement">
      <div>
        <p className="spectral-kicker text-fuchsia-200">Observed cohort state</p>
        <p className="mt-2 max-w-xl text-sm leading-6 text-plum-100">
          Released quality flags are constant in this selected tile. The dashboard preserves the null
          comparison instead of manufacturing a flagged subgroup.
        </p>
      </div>
      <div className="flex gap-7 sm:gap-12">
        <div><p className="font-display text-4xl text-white">{clean.toLocaleString()}</p><p className="mt-1 text-xs uppercase tracking-widest text-fuchsia-200">clean</p></div>
        <div><p className="font-display text-4xl text-white">{flagged.toLocaleString()}</p><p className="mt-1 text-xs uppercase tracking-widest text-fuchsia-200">flagged</p></div>
      </div>
    </div>
  );
}

export default function App() {
  const project = useJson('./project.json');
  const summary = useJson('./results/summary.json');
  const warnings = useJson('./results/warnings.json');
  const benchmarks = useJson('./results/benchmarks.json');

  if (project.loading) {
    return <main className="spectral-page grid min-h-screen place-items-center text-xs uppercase tracking-[0.2em] text-fuchsia-800">Loading spectral audit…</main>;
  }
  if (project.error || !project.data) {
    return <main className="spectral-page grid min-h-screen place-items-center text-red-800">Could not load project.json: {String(project.error)}</main>;
  }

  const p = project.data;
  const metrics = summary.data?.metrics ?? [];
  const isDemo = summary.data?.data_kind === 'synthetic_smoke_test' || summary.data?.data_kind === 'synthetic_demo';

  return (
    <main className="spectral-page min-h-screen">
      <header className="spectral-header">
        <div className="mx-auto max-w-[92rem] px-4 py-4 sm:px-6 lg:px-8">
          <div className="mb-4 flex flex-wrap items-center justify-between gap-3 border-b border-fuchsia-300/20 pb-4 text-xs">
            <div className="flex items-center gap-3 text-fuchsia-100"><Aperture size={17} aria-hidden="true" /><span className="uppercase tracking-[0.2em]">Euclid Q1 · NISP spectroscopy</span></div>
            <div className="flex gap-2">
              <span className="hero-pill">{p.status}</span>
              <span className="hero-pill">Priority {p.priority}/10</span>
            </div>
          </div>
          <div className="grid gap-5 lg:grid-cols-[1.02fr_0.98fr] lg:items-stretch">
            <div className="flex min-h-[30rem] flex-col justify-between py-8 lg:py-14 lg:pr-8">
              <div>
                <p className="spectral-kicker text-fuchsia-300">{p.category}</p>
                <h1 className="font-display mt-5 max-w-4xl text-5xl leading-[0.98] tracking-[-0.045em] text-white sm:text-6xl xl:text-7xl">{p.title}</h1>
                <p className="mt-7 max-w-3xl text-lg leading-8 text-plum-100">{p.question}</p>
              </div>
              <div className="mt-10 flex flex-wrap items-center gap-3 text-xs">
                <span className="rounded-full border border-fuchsia-300/25 bg-fuchsia-950/45 px-4 py-2 text-fuchsia-100">{p.dataMode}</span>
                {summary.data && (
                  <span className={`rounded-full border px-4 py-2 ${isDemo ? 'border-amber-300/40 bg-amber-900/30 text-amber-100' : 'border-emerald-300/35 bg-emerald-900/25 text-emerald-100'}`}>
                    {isDemo ? 'Synthetic demo results' : 'Real data results'}
                  </span>
                )}
              </div>
            </div>
            <LazySpectralHero />
          </div>
        </div>
      </header>

      <div className="mx-auto max-w-[92rem] px-4 pb-14 sm:px-6 lg:px-8">
        {isDemo && (
          <div className="mt-6 flex items-start gap-3 rounded-xl border border-amber-300 bg-amber-50 p-4 text-sm leading-6 text-amber-950">
            <AlertTriangle size={18} className="mt-0.5 shrink-0" aria-hidden="true" />
            These metrics and figures use clearly labelled synthetic demo data, not Euclid observations.
          </div>
        )}

        <CohortStatement metrics={metrics} />

        <section className="mt-14" aria-labelledby="metric-heading">
          <div className="mb-6 flex items-end justify-between gap-4">
            <div><p className="spectral-kicker">Spectral index</p><h2 id="metric-heading" className="font-display mt-1 text-4xl text-plum-950">Recorded measurements</h2></div>
            <p className="hidden max-w-md text-right text-sm leading-6 text-plum-500 md:block">Every reported metric is displayed, including unavailable comparisons from the documented null split.</p>
          </div>
          <div className="metric-grid grid gap-px overflow-hidden rounded-2xl border border-plum-200 bg-plum-200 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4">
            {metrics.map((metric, index) => <MetricCard key={metric.name} metric={metric} index={index} />)}
            {!summary.data && <div className="metric-tile"><p className="font-display text-2xl">No results yet</p><p className="mt-2 text-sm text-plum-500">Run scripts/run_analysis.py first.</p></div>}
          </div>
        </section>

        <section className="mt-14" aria-labelledby="evidence-heading">
          <div className="mb-6"><p className="spectral-kicker">Evidence atlas</p><h2 id="evidence-heading" className="font-display mt-1 text-4xl text-plum-950">Figure gallery</h2></div>
          <div className="evidence-grid">
            {p.figures.map((figure, index) => (
              <figure key={figure.id} className={`evidence-card evidence-card-${index + 1}`}>
                <div className="evidence-image">
                  <img
                    src={`./figures/${figure.id}.svg`}
                    alt={figure.label}
                    className="h-full w-full object-contain"
                    loading={index > 1 ? 'lazy' : 'eager'}
                    onError={(event) => { event.currentTarget.style.display = 'none'; }}
                  />
                </div>
                <figcaption className="flex items-center justify-between gap-3 px-5 py-4">
                  <span className="font-display text-lg text-plum-950">{figure.label}</span>
                  <span className="font-mono text-[0.65rem] text-fuchsia-700">PLATE {String(index + 1).padStart(2, '0')}</span>
                </figcaption>
              </figure>
            ))}
          </div>
        </section>

        <section className="mt-14">
          <Section icon={AlertTriangle} kicker="Live result log" title="Warnings and documented limitations">
            <WarningAudit state={warnings} />
          </Section>
        </section>

        <section className="mt-8 grid gap-8 lg:grid-cols-3">
          <Section icon={ShieldCheck} kicker="Scope control" title="Provenance boundary">
            <p className="text-sm leading-6 text-plum-700">{p.novelty}</p>
            <div className="mt-5 rounded-xl border border-amber-200 bg-amber-50 p-4 text-sm leading-6 text-amber-950">No result is public-ready until validation and provenance checks pass.</div>
            {summary.data?.provenance && (
              <dl className="mt-5 space-y-3 text-xs">
                <div className="flex items-center gap-2"><GitCommit size={14} className="text-fuchsia-700" /><dt className="text-plum-500">git commit</dt><dd className="ml-auto font-mono text-plum-800">{summary.data.provenance.git_commit}</dd></div>
                <div className="flex items-center gap-2"><FileText size={14} className="text-fuchsia-700" /><dt className="text-plum-500">config sha256</dt><dd className="ml-auto max-w-[9rem] truncate font-mono text-plum-800">{summary.data.provenance.config_sha256 ?? 'n/a'}</dd></div>
                <div className="flex items-center gap-2"><Beaker size={14} className="text-fuchsia-700" /><dt className="text-plum-500">package version</dt><dd className="ml-auto font-mono text-plum-800">{summary.data.provenance.package_version}</dd></div>
              </dl>
            )}
          </Section>
          <Section icon={ListChecks} kicker="Acceptance gates" title="Validation contract">
            <ol className="space-y-3 text-sm text-plum-700">
              {p.validationContract.map((item, index) => (
                <li key={item} className="flex gap-3 border-b border-plum-100 pb-3 last:border-0 last:pb-0"><span className="font-mono text-fuchsia-700">{String(index + 1).padStart(2, '0')}</span><span>{item}</span></li>
              ))}
            </ol>
          </Section>
          <ConfidenceExplorer metrics={metrics} />
        </section>

        <section className="method-band mt-14 grid gap-8 lg:grid-cols-[1.15fr_0.85fr]">
          <Section icon={Sparkles} kicker="Analysis chain" title="Methodology" dark>
            <p className="text-sm leading-7 text-plum-100">{p.methodology}</p>
          </Section>
          <Section icon={BookOpen} kicker="Interpretive boundary" title="Assumptions and limitations" dark>
            <p className="spectral-kicker text-fuchsia-300">Assumptions</p>
            <ul className="mt-3 space-y-3 text-sm leading-6 text-plum-100">
              {p.assumptions.map((assumption) => <li key={assumption} className="border-l border-fuchsia-400/50 pl-3">{assumption}</li>)}
            </ul>
            <p className="spectral-kicker mt-7 text-amber-300">Limitations</p>
            <ul className="mt-3 space-y-3 text-sm leading-6 text-plum-100">
              {p.limitations.map((limitation) => <li key={limitation} className="border-l border-amber-400/50 pl-3">{limitation}</li>)}
            </ul>
          </Section>
        </section>

        <footer className="mt-8 grid gap-8 lg:grid-cols-[1.2fr_0.8fr]">
          <Section icon={Download} kicker="Reproducibility" title="Downloads and provenance manifest">
            <div className="flex flex-wrap gap-2 text-sm">
              <a className="download-link" href="./manifest.csv" download>data/manifest.csv</a>
              <a className="download-link" href="./results/summary.json" download>results/summary.json</a>
              {benchmarks.data && <a className="download-link" href="./results/benchmarks.json" download>results/benchmarks.json</a>}
            </div>
            <p className="mt-5 text-xs leading-5 text-plum-500">The manifest records product identifiers, archive source, retrieval time, checksums, file sizes, selection reasons and usage terms.</p>
          </Section>
          <Section icon={Database} kicker="Credit" title="Citation and licence">
            <p className="text-sm text-plum-700">Author: {p.citation.author}</p>
            <p className="mt-2 text-sm text-plum-700">Licence: {p.citation.license}</p>
            <a className="mt-4 inline-block text-sm text-fuchsia-800 underline-offset-4 hover:underline" href={p.citation.repository}>{p.citation.repository}</a>
          </Section>
        </footer>
      </div>
    </main>
  );
}
