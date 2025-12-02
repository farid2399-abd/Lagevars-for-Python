# estimation_largevar.py
import numpy as np
from scipy.linalg import eig, pinv
from typing import Tuple, Optional
import warnings

def estimar_parametros_largevar(
    X: np.ndarray,
    r: int = 1,
    normalize: bool = True,
    demean: bool = True
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Estima parámetros beta y alpha usando el método Largevar
    
    Parámetros:
    -----------
    X : np.ndarray
        Datos de forma (T, N)
    r : int
        Rango de cointegración a estimar
    normalize : bool
        Si True, normaliza beta para que el primer elemento sea 1
    demean : bool
        Si True, elimina la media de los datos
    
    Retorna:
    --------
    beta_est : np.ndarray o None
        Vector beta estimado, forma (N, r)
    alpha_est : np.ndarray o None
        Matriz alpha estimada, forma (N, r)
    
    Ejemplo:
    --------
    >>> X, _, _, _ = simular_vecm_configurable(T=100, N=5)
    >>> beta, alpha = estimar_parametros_largevar(X, r=1)
    """
    try:
        T, N = X.shape
        
        if T < 2:
            raise ValueError("Se necesitan al menos 2 observaciones")
        
        t = T - 1
        
        # Paso 1: Transformación Largevar (eliminar tendencia estocástica)
        X_tilde = np.zeros((N, t))
        dX = np.zeros((N, t))
        
        for i in range(t):
            # Transformación de ecuaciones (4) y (5)
            X_tilde[:, i] = X[i, :] - (i / t) * (X[-1, :] - X[0, :])
            dX[:, i] = X[i + 1, :] - X[i, :]
        
        # Paso 2: Centrar datos si es necesario
        if demean:
            R0 = (dX - dX.mean(axis=1, keepdims=True)).T
            R1 = (X_tilde - X_tilde.mean(axis=1, keepdims=True)).T
        else:
            R0 = dX.T
            R1 = X_tilde.T
        
        # Paso 3: Matrices de covarianza
        S00 = (R0.T @ R0) / t
        Skk = (R1.T @ R1) / t
        S0k = (R0.T @ R1) / t
        
        # Paso 4: Problema de valores propios generalizado
        # Resolver: S0k * S00^{-1} * S0k' * v = λ * Skk * v
        inv_S00 = pinv(S00)
        M = pinv(Skk) @ S0k @ inv_S00 @ S0k.T
        
        # Valores propios y vectores
        eigenvals, eigenvecs = eig(M)
        
        # Ordenar por valores propios descendentes
        idx = np.argsort(np.real(eigenvals))[::-1]
        eigenvals_sorted = np.real(eigenvals[idx])
        eigenvecs_sorted = np.real(eigenvecs[:, idx])
        
        # Seleccionar r vectores propios para beta
        beta_est = eigenvecs_sorted[:, :r]
        
        # Paso 5: Estimar alpha por MCO
        if r > 0:
            Z = R1 @ beta_est
            alpha_est = pinv(Z.T @ Z) @ Z.T @ R0
            alpha_est = alpha_est.T
        else:
            alpha_est = np.zeros((N, r))
        
        # Paso 6: Normalización
        if normalize and r > 0 and beta_est.shape[0] > 0:
            for i in range(r):
                if beta_est[0, i] != 0:
                    factor = beta_est[0, i]
                    beta_est[:, i] = beta_est[:, i] / factor
                    alpha_est[:, i] = alpha_est[:, i] * factor
        
        # Información adicional
        info = {
            'valores_propios': eigenvals_sorted,
            'varianza_explicada': eigenvals_sorted[:r] / np.sum(eigenvals_sorted),
            'T': T,
            'N': N,
            'r': r
        }
        
        return beta_est, alpha_est, info
        
    except Exception as e:
        warnings.warn(f"Error en estimación Largevar: {str(e)}")
        return None, None, {}


def comparar_metodos_estimacion(
    X: np.ndarray,
    r: int = 1,
    p: int = 0
) -> dict:
    """
    Compara métodos Johansen y Largevar
    
    Parámetros:
    -----------
    X : np.ndarray
        Datos de entrada
    r : int
        Rango de cointegración
    p : int
        Orden del VECM
    
    Retorna:
    --------
    comparacion : dict
        Diccionario con resultados comparativos
    """
    # Estimación Johansen
    beta_j, alpha_j = estimar_parametros_johansen(X, r=r, p=p)
    
    # Estimación Largevar
    beta_l, alpha_l, info_l = estimar_parametros_largevar(X, r=r)
    
    comparacion = {
        'johansen': {
            'beta': beta_j,
            'alpha': alpha_j,
            'metodo': 'Johansen MLE'
        },
        'largevar': {
            'beta': beta_l,
            'alpha': alpha_l,
            'info': info_l,
            'metodo': 'Largevar Eigenvalue'
        }
    }
    
    # Calcular similitud entre betas si ambos están disponibles
    if beta_j is not None and beta_l is not None:
        # Correlación entre betas
        if beta_j.ndim == 2 and beta_l.ndim == 2:
            beta_j_flat = beta_j[:, 0] if beta_j.shape[1] > 0 else beta_j.flatten()
            beta_l_flat = beta_l[:, 0] if beta_l.shape[1] > 0 else beta_l.flatten()
            
            correlacion = np.corrcoef(beta_j_flat, beta_l_flat)[0, 1]
            distancia = np.linalg.norm(beta_j_flat - beta_l_flat)
            
            comparacion['similitud'] = {
                'correlacion': correlacion,
                'distancia_euclidiana': distancia
            }
    
    return comparacion
