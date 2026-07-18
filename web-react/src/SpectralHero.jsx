import { useEffect, useRef, useState } from 'react';
import { Canvas, useFrame, useThree } from '@react-three/fiber';

function useReducedMotion() {
  const [reduced, setReduced] = useState(false);

  useEffect(() => {
    const query = window.matchMedia('(prefers-reduced-motion: reduce)');
    const update = () => setReduced(query.matches);
    update();
    query.addEventListener('change', update);
    return () => query.removeEventListener('change', update);
  }, []);

  return reduced;
}

function FrameTicker({ active }) {
  const invalidate = useThree((state) => state.invalidate);

  useEffect(() => {
    invalidate();
    if (!active) return undefined;
    const timer = window.setInterval(invalidate, 80);
    return () => window.clearInterval(timer);
  }, [active, invalidate]);

  return null;
}

const beamColors = ['#f9a8d4', '#e879f9', '#c084fc', '#67e8f9', '#fde68a'];

function DispersedBeam() {
  return (
    <group position={[1.2, 0.15, 0]} rotation={[0, 0, -0.12]}>
      {beamColors.map((color, index) => (
        <mesh key={color} position={[0.55 + index * 0.15, (index - 2) * 0.09, 0]} rotation={[0, 0, (index - 2) * 0.022]}>
          <boxGeometry args={[1.8, 0.045, 0.035]} />
          <meshBasicMaterial color={color} transparent opacity={0.44 - index * 0.025} />
        </mesh>
      ))}
    </group>
  );
}

function EuclidNispModel({ animate }) {
  const model = useRef(null);

  useFrame((state, delta) => {
    if (!model.current || !animate) return;
    model.current.rotation.y += delta * 0.14;
    model.current.position.y = Math.sin(state.clock.elapsedTime * 0.5) * 0.07;
  });

  return (
    <group ref={model} rotation={[0.2, -0.65, -0.08]} position={[-0.45, 0, 0]}>
      <mesh>
        <cylinderGeometry args={[0.74, 0.9, 1.85, 16]} />
        <meshStandardMaterial color="#ddd6e6" roughness={0.45} metalness={0.55} />
      </mesh>
      <mesh position={[0, 1, 0]}>
        <cylinderGeometry args={[0.73, 0.73, 0.12, 24]} />
        <meshStandardMaterial color="#4a164f" roughness={0.35} metalness={0.62} />
      </mesh>
      <mesh position={[0, 1.07, 0]} rotation={[Math.PI / 2, 0, 0]}>
        <circleGeometry args={[0.57, 32]} />
        <meshStandardMaterial color="#13091b" roughness={0.18} metalness={0.45} />
      </mesh>
      <mesh position={[0, -1.08, 0]}>
        <boxGeometry args={[1.28, 0.38, 1.15]} />
        <meshStandardMaterial color="#6b446e" roughness={0.42} metalness={0.55} />
      </mesh>
      <mesh position={[0, -1.32, 0]} rotation={[0.05, 0, 0]}>
        <cylinderGeometry args={[1.9, 1.9, 0.08, 6]} />
        <meshStandardMaterial color="#b97886" roughness={0.58} metalness={0.42} side={2} />
      </mesh>
      <mesh position={[0, -1.25, 0]}>
        <cylinderGeometry args={[0.48, 0.48, 0.16, 12]} />
        <meshStandardMaterial color="#29122f" roughness={0.4} metalness={0.6} />
      </mesh>
      <mesh position={[0.62, -0.65, 0.28]} rotation={[0, 0.3, 0]}>
        <boxGeometry args={[0.48, 0.52, 0.18]} />
        <meshStandardMaterial color="#9d174d" emissive="#701a75" emissiveIntensity={0.24} roughness={0.38} metalness={0.5} />
      </mesh>
      <DispersedBeam />
    </group>
  );
}

export default function SpectralHero() {
  const reducedMotion = useReducedMotion();

  return (
    <figure className="relative h-full w-full" aria-label="Procedural Euclid-like NISP spectral instrument illustration">
      <Canvas
        camera={{ position: [5.4, 3.2, 6.7], fov: 38 }}
        dpr={[1, 1.5]}
        frameloop="demand"
        gl={{ antialias: false, alpha: true, powerPreference: 'low-power' }}
      >
        <FrameTicker active={!reducedMotion} />
        <ambientLight intensity={1.2} />
        <directionalLight position={[5, 6, 4]} intensity={2.4} color="#fdf4ff" />
        <pointLight position={[-3, -1, 4]} intensity={1.2} color="#f0abfc" />
        <EuclidNispModel animate={!reducedMotion} />
      </Canvas>
      <div className="pointer-events-none absolute inset-x-0 bottom-0 flex items-end justify-between gap-4 bg-gradient-to-t from-[#1b081e] via-[#1b081e]/85 to-transparent px-5 pb-5 pt-16">
        <div>
          <p className="spectral-kicker text-fuchsia-300">NISP slitless spectral path</p>
          <p className="mt-1 text-xs text-plum-200">Procedural low-poly mission schematic</p>
        </div>
        <figcaption className="max-w-40 text-right text-[0.65rem] uppercase tracking-[0.12em] text-fuchsia-100/80">
          Stylized illustration, not flight data
        </figcaption>
      </div>
    </figure>
  );
}
