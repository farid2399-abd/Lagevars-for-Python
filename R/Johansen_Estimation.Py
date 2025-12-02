# estimation_johansen.py
import numpy as np
from statsmodels.tsa.vector_ar.vecm import VECM
from typing import Tuple, Optional
import warnings

def estimar_parametros_johansen(
    X: np.ndarray,
    r: int = 1,
    p: int = 0,
    normalize: bool = True,
    deterministic: str = 'n'
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Estima parámetros beta y alpha usando el método de Johansen (MLE)
    
    Parámetros:
    -----------
    X : np.ndarray
        Datos de forma (T, N)
    r : int
        Rango de cointegración a estimar
    p : int
        Orden del VAR en niveles (para VECM, p-1 retardos en diferencias)
    normalize : bool
        Si True, normaliza beta para que el primer elemento sea 1
    deterministic : str
        Componente determinístico ('n': ninguno, 'co': constante en cointegración)
    
    Retorna:
    --------
    beta_est : np.ndarray o None
        Vector beta estimado, forma (N, r)
    alpha_est : np.ndarray o None
        Matriz alpha estimada, forma (N, r)
    
    Ejemplo:
    --------
    >>> X, _, _, _ = simular_vecm_configurable(T=100, N=5)
    >>> beta, alpha = estimar_parametros_johansen(X, r=1, p=0)
    """
    try:
        T, N = X.shape
        
        # Verificar dimensiones
        if r > N:
            raise ValueError(f"El rango r={r} no puede ser mayor que N={N}")
        
        if T < N:
            warnings.warn(f"Pocas observaciones (T={T}) para N={N} series")
        
        # Ajustar k_ar_diff para statsmodels (número de retardos en diferencias)
        k_ar_diff = max(p - 1, 0) if p > 0 else 0
        
        # Estimar con VECM de statsmodels
        model = VECM(
            X,
            k_ar_diff=k_ar_diff,
            coint_rank=r,
            deterministic=deterministic
        )
        result = model.fit()
        
        beta_est = result.beta
        alpha_est = result.alpha
        
        # Normalización
        if normalize and beta_est is not None:
            for i in range(min(r, beta_est.shape[1])):
                if beta_est[0, i] != 0:
                    factor = beta_est[0, i]
                    beta_est[:, i] = beta_est[:, i] / factor
                    alpha_est[:, i] = alpha_est[:, i] * factor
        
        return beta_est, alpha_est
        
    except Exception as e:
        warnings.warn(f"Error en estimación Johansen: {str(e)}")
        return None, None


def comparar_estimacion_johansen(
    X: np.ndarray,
    rangos: list = [1, 2, 3],
    p: int = 0
) -> dict:
    """
    Compara estimaciones Johansen para diferentes rangos
    
    Parámetros:
    -----------
    X : np.ndarray
        Datos de entrada
    rangos : list
        Lista de rangos a probar
    p : int
        Ordel del VECM
    
    Retorna:
    --------
    resultados : dict
        Diccionario con resultados por rango
    """
    resultados = {}
    
    for r in rangos:
        beta, alpha = estimar_parametros_johansen(X, r=r, p=p)
        
        if beta is not None:
            resultados[r] = {
                'beta': beta,
                'alpha': alpha,
                'rango': r,
                'beta_norm': np.linalg.norm(beta) if beta is not None else None,
                'alpha_norm': np.linalg.norm(alpha) if alpha is not None else None
            }
    
    return resultados
